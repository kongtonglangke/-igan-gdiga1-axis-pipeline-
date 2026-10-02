# -*- coding: utf-8 -*-
"""主证据2 补数据: 抓取 chr1 (GALNT2) 与 chr17 (ST6GALNAC2) 的 region_cache
复用主证据1 stepA 已验证的远程 BGZF/tabix 抓取管线 (5 QTD 语境)。
Ensembl 坐标 (GRCh38):
  GALNT2    ENSG00000143641  chr1 : 230,057,990-230,282,128
  ST6GALNAC2 ENSG00000070731 chr17: 76,564,424-76,586,956
窗口 = 基因 start/end 各外扩 500kb (cis-eQTL 常规窗口)。
"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import os, sys, time, urllib.request, concurrent.futures as cf

sys.path.insert(0, _paths.LIB)   # 随包分发的 tabix_remote.py
from tabix_remote import parse_tbi, region_file_ranges, decompress_bgzf

TBI_DIR = _paths.raw("eQTL_Catalogue", "tbi")   # .tbi 索引缓存（缺失时自动下载）
CACHE_DIR = _paths.raw("eQTL_Catalogue", "region_cache")
os.makedirs(CACHE_DIR, exist_ok=True)

def tbi_path(qtd, qts):
    """{qtd}.tbi 索引路径；缺失时自动从 eQTL Catalogue 下载（第三方索引，不随仓库分发）。"""
    os.makedirs(TBI_DIR, exist_ok=True)
    p = os.path.join(TBI_DIR, f"{qtd}.tbi")
    if not os.path.exists(p):
        url = f"{BASE}/{qts}/{qtd}/{qtd}.all.tsv.gz.tbi"
        print(f"[tbi] 下载索引 {url}", flush=True)
        urllib.request.urlretrieve(url, p)
    return p


BASE = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats"
DATASETS = [
    ("OneK1K_B_naive",        "QTD000608", "QTS000038"),
    ("OneK1K_B_memory",       "QTD000607", "QTS000038"),
    ("OneK1K_B_intermediate", "QTD000606", "QTS000038"),
    ("GTEx_v10_blood",        "QTD000356", "QTS000015"),
    ("GTEx_v10_kidney_cortex","QTD000261", "QTS000015"),
]
# (chrom, beg, end) 1-based 闭区间
REGIONS = [
    ("1",  229_557_990, 230_782_128),  # GALNT2 ±500kb
    ("17", 76_064_424,  77_086_956),   # ST6GALNAC2 ±500kb
]


def file_size(url):
    req = urllib.request.Request(url, headers={"Range": "bytes=0-0", "User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=120) as r:
        cr = r.headers.get("Content-Range", "")
        return int(cr.split("/")[1]) if "/" in cr else None


def fetch_blob(url, fbeg, fend_excl, total=None):
    if total is not None and fend_excl > total:
        fend_excl = total
    req = urllib.request.Request(url, headers={"Range": f"bytes={fbeg}-{fend_excl-1}", "User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=300) as r:
        if r.status == 200 and fbeg != 0:
            raise RuntimeError("server ignored Range")
        return r.read()


def fetch_region_rows(url, tbi, chrom, beg, end):
    found, ranges = region_file_ranges(tbi, chrom, beg, end)
    if not found:
        return []
    total = file_size(url)
    parts = []
    for fbeg, fend in ranges:
        raw = fetch_blob(url, fbeg, fend + 65536, total=total)
        parts.append(decompress_bgzf(raw))
    return b"".join(parts).split(b"\n")


def cache_one(qtd, qts, chrom, beg, end):
    out = os.path.join(CACHE_DIR, f"{qtd}.{chrom}.tsv")
    if os.path.exists(out) and os.path.getsize(out) > 1000:
        print(f"[skip] {out} exists", flush=True)
        return
    tbi = parse_tbi(tbi_path(qtd, qts))
    url = f"{BASE}/{qts}/{qtd}/{qtd}.all.tsv.gz"
    t0 = time.time()
    rows = fetch_region_rows(url, tbi, chrom, beg, end)
    with open(out, "wb") as f:
        for r in rows:
            if r.strip():
                f.write(r + b"\n")
    print(f"[OK] {qtd} {chrom}:{beg}-{end} rows={len(rows)} {time.time()-t0:.0f}s -> {out}", flush=True)


def job(qtd_qts):
    qtd, qts = qtd_qts
    for chrom, beg, end in REGIONS:
        cache_one(qtd, qts, chrom, beg, end)


if __name__ == "__main__":
    t_start = time.time()
    with cf.ThreadPoolExecutor(max_workers=4) as ex:
        futs = [ex.submit(job, (q, qs)) for _, q, qs in DATASETS]
        for f in cf.as_completed(futs):
            f.result()
    print(f"ALL DONE in {time.time()-t_start:.0f}s")
    for fn in sorted(os.listdir(CACHE_DIR)):
        if fn.endswith((".1.tsv", ".17.tsv")):
            print("  ", fn, os.path.getsize(os.path.join(CACHE_DIR, fn)))
