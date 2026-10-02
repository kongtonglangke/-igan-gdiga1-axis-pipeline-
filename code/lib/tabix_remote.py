# -*- coding: utf-8 -*-
"""纯 Python 远程 BGZF/tabix 区域读取器 v2（pysam 的 Windows 替代，实测修正版）

背景：eQTL Catalogue .all.tsv.gz（每数据集 1.8-2GB）含全 nominal cis-eQTL 关联，
主证据1 需要"lead SNP 本尊 × 细胞型"状态，只能查 all 文件。沙箱无 pysam/编译器，
且全文件下载不划算 → 解析 .tbi 索引 + HTTP Range 只拉目标区域所在 BGZF 块。

.v2 相对 v1 的修正（2026-09-02 实测 .tbi 头字节）：
  .tbi 头 = magic(4) + 8×int32: n_ref=23, fmt=0, col_seq=2, col_beg=3, col_end=3,
            meta=35('#'), skip=1, l_nm=59，随后 l_nm 字节 NUL 分隔的染色体名（23 个，含 X）。
  旧 v1 误把第 5 字节当 1 字节 fmt → 整体错位 → UnicodeDecodeError。

BGZF: 每块为独立 gzip member（压缩后≤64KB）；tabix 虚拟偏移 voffset =
(块文件偏移 << 16) | 块内偏移。.tbi 的 bin 索引给出 region -> chunk(voffset对)。
坐标排序的单点型数据（variant position）只落 16kb 叶 bin，祖先 bin 一般 n_chunk=0，
因此 region 查询的 chunk 天然局部化，不会误抓整条染色体。

本模块只依赖标准库：struct / urllib.request / gzip / io。
"""
import struct, urllib.request, gzip, io, sys

CHROM_NORM = {}
for i in range(1, 23):
    CHROM_NORM[str(i)] = str(i)
CHROM_NORM.update({"X": "X", "23": "X", "chrX": "X", "chr7": "7", "chr9": "9", "chr1": "1", "chr17": "17"})

_HTTP = {"User-Agent": "Mozilla/5.0 (research-eqtl-tabix-reader)"}


def _i32(d, o):
    return struct.unpack_from("<i", d, o)[0]


def parse_tbi(path):
    """解析 .tbi → dict(refs 有序 list, per-ref bins 与 linear)"""
    with open(path, "rb") as f:
        d = f.read()
    assert d[:4] == b"TBI\x01", "not a tabix index"
    o = 4
    n_ref = _i32(d, o); o += 4
    fmt = _i32(d, o); o += 4
    col_seq = _i32(d, o); o += 4
    col_beg = _i32(d, o); o += 4
    col_end = _i32(d, o); o += 4
    meta = _i32(d, o); o += 4
    skip = _i32(d, o); o += 4
    l_nm = _i32(d, o); o += 4
    names = d[o:o + l_nm].rstrip(b"\x00").split(b"\x00")
    names = [x.decode("utf-8") for x in names if x]
    o += l_nm
    refs = []
    for ri in range(n_ref):
        n_bin = _i32(d, o); o += 4
        bins = {}
        for _bi in range(n_bin):
            bid = struct.unpack_from("<I", d, o)[0]; o += 4
            n_chunk = _i32(d, o); o += 4
            chunks = []
            for _ci in range(n_chunk):
                vb, ve = struct.unpack_from("<QQ", d, o); o += 16
                chunks.append((vb, ve))
            bins[bid] = chunks
        n_intv = _i32(d, o); o += 4
        intv = list(struct.unpack_from("<%dQ" % n_intv, d, o)) if n_intv else []
        o += 8 * n_intv
        refs.append({"name": names[ri], "bins": bins, "intv": intv})
    # 索引尾部：n_no_coor (uint64, 通常=0)，消费掉使解析吃满整个文件
    if o + 8 <= len(d):
        n_no_coor = struct.unpack_from("<Q", d, o)[0]; o += 8
    else:
        n_no_coor = None
    # 校验：若 o 恰好==len(d) 说明头与索引字段长度全部猜对
    consumed = o
    return {"n_ref": n_ref, "fmt": fmt, "col_seq": col_seq, "col_beg": col_beg,
            "col_end": col_end, "meta": meta, "skip": skip, "names": names,
            "refs": refs, "header_bytes": consumed, "file_bytes": len(d)}


def reg2bins(beg, end):
    """UCSC reg2bins：返回覆盖 [beg,end)（0-based half-open）的全部候选 bin id"""
    lst = [0]
    end -= 1
    k = 1 + (beg >> 26)
    while k <= 1 + (end >> 26):
        lst.append(k); k += 1
    k = 9 + (beg >> 23)
    while k <= 9 + (end >> 23):
        lst.append(k); k += 1
    k = 73 + (beg >> 20)
    while k <= 73 + (end >> 20):
        lst.append(k); k += 1
    k = 585 + (beg >> 17)
    while k <= 585 + (end >> 17):
        lst.append(k); k += 1
    k = 4681 + (beg >> 14)
    while k <= 4681 + (end >> 14):
        lst.append(k); k += 1
    return sorted(set(lst))


def region_file_ranges(tbi, chrom, beg, end):
    """给定 1-based 区间 [beg,end]，返回需要抓取的 BGZF 文件字节区间列表
    返回 (ref_found: bool, [(fbeg, fend_excl), ...])，fbeg 均为块边界。
    策略：取覆盖区间的候选 bin（叶+祖先）中有数据的 bin 的 chunk voffset，
    全部换算成块偏移后按连续性合并（相邻/重叠合并），宁可少跨块也不全抓。
    """
    idx = None
    for i, r in enumerate(tbi["refs"]):
        if r["name"] == chrom or CHROM_NORM.get(r["name"]) == chrom:
            idx = i; break
    if idx is None:
        return False, []
    ref = tbi["refs"][idx]
    beg0, end0 = beg - 1, end  # 转 0-based half-open
    cand = reg2bins(beg0, end0)
    spans = []  # 文件块区间 (块号), 用块号(>>16)
    for bid in cand:
        for vb, ve in ref["bins"].get(bid, []):
            spans.append((vb >> 16, ve >> 16))
    if not spans:
        return True, []
    spans.sort()
    merged = []
    cur = list(spans[0])
    for s in spans[1:]:
        if s[0] <= cur[1]:          # 相邻或重叠块号 → 合并
            cur[1] = max(cur[1], s[1])
        else:
            merged.append(tuple(cur)); cur = list(s)
    merged.append(tuple(cur))
    # 块号 → 字节偏移（BGZF 每块在文件中是连续字节，块号即文件偏移）
    return True, [(b * 1, e * 1) for b, e in merged]


def _file_size(url):
    """用 bytes=0-0 探文件总长（HEAD 被代理吞掉时改用小 Range GET）"""
    req = urllib.request.Request(url, headers={**{"Range": "bytes=0-0"}, **_HTTP})
    with urllib.request.urlopen(req, timeout=120) as r:
        cr = r.headers.get("Content-Range", "")
        if "/" in cr:
            return int(cr.split("/")[1])
        return None


def fetch_bytes(url, fbeg, fend_excl, total=None):
    """Range 抓取 [fbeg, fend_excl)。若服务器忽略 Range(200) 则熔断，避免下整文件。"""
    if total is not None and fend_excl > total:
        fend_excl = total
    req = urllib.request.Request(url, headers={**{"Range": f"bytes={fbeg}-{fend_excl - 1}"}, **_HTTP})
    with urllib.request.urlopen(req, timeout=180) as r:
        if r.status == 200 and fbeg != 0:
            raise RuntimeError(f"server ignored Range (status 200) for {url}")
        return r.read()


def decompress_bgzf(data):
    """把连续 BGZF 块解压为字节流（每块独立 gzip member，按 BSIZE 逐块解压拼接）"""
    out = io.BytesIO()
    pos, n = 0, len(data)
    while pos + 18 <= n:
        if data[pos] != 0x1F or data[pos + 1] != 0x8B:
            break
        bsize = struct.unpack_from("<H", data, pos + 16)[0] + 1  # BSIZE = 总块长-1
        block = data[pos:pos + bsize]
        if len(block) < bsize:
            break
        try:
            out.write(gzip.decompress(block))
        except Exception:
            break
        pos += bsize
    return out.getvalue()


def fetch_region_lines(url, tbi, chrom, beg, end, max_rows=10_000_000):
    """远程抓取区域原始行（bytes 行列表），供调用方自行 split 解析。
    beg/end 为 1-based 闭区间坐标（仅用于过滤，实际过滤交给调用方）。"""
    found, ranges = region_file_ranges(tbi, chrom, beg, end)
    if not found:
        return []
    total = _file_size(url)
    parts = []
    for fbeg, fend in ranges:
        raw = fetch_bytes(url, fbeg, fend + 65536, total=total)  # 尾部+1块缓冲
        parts.append(decompress_bgzf(raw))
    blob = b"".join(parts)
    return blob.split(b"\n")


if __name__ == "__main__":
    # 用法: python tabix_remote.py <tbi> <chrom> <beg> <end> <url>
    # 例: python tabix_remote.py QTD000608.tbi 7 6707215 7707215 <url>
    tbi_path = sys.argv[1]
    chrom = sys.argv[2]
    beg, end = int(sys.argv[3]), int(sys.argv[4])
    url = sys.argv[5]
    tbi = parse_tbi(tbi_path)
    print(f"tbi header consumed {tbi['header_bytes']}/{tbi['file_bytes']} bytes, "
          f"n_ref={tbi['n_ref']} names({len(tbi['names'])})={tbi['names'][:6]}...{tbi['names'][-3:]}")
    print(f"query {chrom}:{beg}-{end}")
    found, ranges = region_file_ranges(tbi, chrom, beg, end)
    print(f"ref found={found}, file byte ranges={ranges}")
    if found and ranges:
        lines = fetch_region_lines(url, tbi, chrom, beg, end)
        print(f"total lines fetched (incl header/blanks): {len(lines)}")
        for ln in lines[:6]:
            print("  |", ln.decode("utf-8", "replace")[:200])
