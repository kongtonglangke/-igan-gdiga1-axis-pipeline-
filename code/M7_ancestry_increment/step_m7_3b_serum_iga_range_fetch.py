# -*- coding: utf-8 -*-
"""远端提取 chr7(C1GALT1)/chr9(GALNT12) 基因区（hg19 ±500 kb）。

关键事实（2026-09-26 实证）：META.IGA.LEVELS 文件的排序键是 **SNP 字符串 `CHR:BP`** 的
字典序（非数值 BP）。因此定位必须按 SNP 字符串二分；字节跨度 [lo_str, hi_str] 内的
记录按数值 BP 过滤即可（7 位 BP 的字符串必然落在这两个端点之间）。
"""
import ssl, time, urllib.request, os, sys
from concurrent.futures import ThreadPoolExecutor

CTX = ssl.create_default_context(); CTX.check_hostname = False; CTX.verify_mode = ssl.CERT_NONE
UA = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'}
URL = 'https://www.columbiamedicine.org/divisions/kiryluk/gwas/IgA/META.IGA.LEVELS.ALL.COMBINED.txt'
TOTAL = 596831842
OUT = 'data/igaLevels_remote_part.tsv'
TARGETS = {'C1GALT1': ('7', 7196565, 7288282), 'GALNT12': ('9', 101569981, 101612363)}
PAD = 500000

def log(*a):
    print('[%s]' % time.strftime('%H:%M:%S'), *a, flush=True)

def fetch(a, n, tries=4, timeout=25):
    last = None
    for t in range(tries):
        try:
            h = dict(UA); h['Range'] = 'bytes=%d-%d' % (a, min(a + n - 1, TOTAL - 1))
            with urllib.request.urlopen(urllib.request.Request(URL, headers=h),
                                        timeout=timeout, context=CTX) as f:
                d = f.read()
            if len(d) >= 16:
                return d
            last = 'short %d' % len(d)
        except Exception as e:
            last = type(e).__name__
            time.sleep(0.3 + 0.25 * t)
    raise RuntimeError('fetch %d+%d: %s' % (a, n, last))

def fetch_par(a, n, chunk=262144, workers=8):
    nch = (n + chunk - 1) // chunk
    out = [None] * nch
    def job(i):
        out[i] = fetch(a + i * chunk, min(chunk, n - i * chunk))
    with ThreadPoolExecutor(workers) as ex:
        list(ex.map(job, range(nch)))
    return b''.join(out)

def first_snp(off, n=8192):
    """窗口内首个「完整行」的 SNP 字符串（跳过被截断的首行）"""
    d = fetch(off, n)
    lines = d.split(b'\n')
    if off != 0:
        lines = lines[1:]
    for ln in lines:
        p = ln.split(b'\t')
        if len(p) >= 3 and len(p[0]) >= 3:
            return p[0].decode()
    return None

def locate_str(snp):
    """首个 SNP 字符串 >= snp 的字节偏移（文件按 SNP 字符串字典序排）"""
    lo, hi = 0, TOTAL - 1
    steps = 0
    while lo < hi:
        mid = (lo + hi) // 2
        s = None
        for _ in range(4):
            try:
                s = first_snp(mid); break
            except Exception as e:
                log('    probe retry @%d %s' % (mid, e)); time.sleep(0.8)
        if s is None:
            raise RuntimeError('probe fail @%d' % mid)
        steps += 1
        if s >= snp:
            hi = mid
        else:
            lo = mid + 1
    log('    locate "%s" -> %d (%d steps)' % (snp, lo, steps))
    return lo

def grab(g, chrom, s, e):
    bp_lo, bp_hi = s - PAD, e + PAD
    t0 = time.time()
    lo_off = max(0, locate_str('%s:%d' % (chrom, bp_lo)) - 1)
    hi_off = locate_str('%s:%d' % (chrom, bp_hi))
    n = min(TOTAL - lo_off, hi_off - lo_off + 65536)
    log('  %s 取 [%d..%d] %d B (%.2f MB)' % (g, lo_off, hi_off, n, n / 1048576.0))
    buf = fetch_par(lo_off, n)
    if lo_off > 0:
        i = buf.find(b'\n')
        buf = buf[i + 1:] if i >= 0 else b''
    i = buf.rfind(b'\n')
    if i >= 0:
        buf = buf[:i + 1]
    keep, bps = [], []
    for ln in buf.split(b'\n'):
        p = ln.split(b'\t')
        if len(p) >= 8 and p[1].decode().strip() == chrom:
            b = int(p[2])
            if bp_lo <= b <= bp_hi:
                keep.append(ln); bps.append(b)
    log('  %s 完成 %d 条, BP %d..%d, %.1fs' %
        (g, len(keep), min(bps) if bps else -1, max(bps) if bps else -1, time.time() - t0))
    return g, b''.join(x + b'\n' for x in keep), (min(bps), max(bps), len(keep))

def main():
    parts, summary = {}, {}
    def run(g):
        c, s, e = TARGETS[g]
        _, blob, info = grab(g, c, s, e)
        parts[g] = blob; summary[g] = info
    with ThreadPoolExecutor(2) as ex:
        list(ex.map(run, list(TARGETS)))
    with open(OUT, 'wb') as f:
        f.write(b'SNP\tCHR\tBP_hg19\tA1\tA2\tBETA\tSE\tP\n')
        for g in ('C1GALT1', 'GALNT12'):
            f.write(parts[g])
    for g in ('C1GALT1', 'GALNT12'):
        lo, hi, n = summary[g]
        want_lo = TARGETS[g][1] - PAD; want_hi = TARGETS[g][2] + PAD
        ok = (n > 0) and (lo <= want_lo + 200000) and (hi >= want_hi - 200000)
        log('  %s 覆盖 BP %d..%d（目标 %d..%d）%s' % (g, lo, hi, want_lo, want_hi, 'OK' if ok else '⚠ 需检查'))
    log('saved -> %s %d bytes' % (OUT, os.path.getsize(OUT)))

if __name__ == '__main__':
    main()
