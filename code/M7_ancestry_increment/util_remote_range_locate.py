# -*- coding: utf-8 -*-
"""远端 Range 定位器：对 Kiryluk 汇总统计（按 CHR 分组、组内按 BP 升序）建粗索引，
再二分定位目标染色体区间内的坐标 → 只下载所需区段。
"""
import ssl, time, urllib.request, sys, json, os
from concurrent.futures import ThreadPoolExecutor

CTX = ssl.create_default_context(); CTX.check_hostname = False; CTX.verify_mode = ssl.CERT_NONE
UA = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'}

def fetch_range(url, a, n, tries=8, timeout=90):
    need = n
    for t in range(tries):
        try:
            h = dict(UA); h['Range'] = 'bytes=%d-%d' % (a, a + n - 1)
            with urllib.request.urlopen(urllib.request.Request(url, headers=h),
                                        timeout=timeout, context=CTX) as f:
                d = f.read()
            if len(d) >= min(n, 16):
                return d
        except Exception:
            time.sleep(0.5 + 0.4 * t)
    raise RuntimeError('range fetch failed at %d' % a)

def probe(url, off, n=16384):
    """返回该窗口内【完整行】的 (chr, bp) 列表"""
    d = fetch_range(url, off, n)
    lines = d.split(b'\n')
    if off != 0:
        lines = lines[1:]           # 首行可能被截断
    if not d.endswith(b'\n'):
        lines = lines[:-1]          # 末行可能被截断
    out = []
    for ln in lines:
        p = ln.split(b'\t')
        if len(p) >= 3:
            try:
                out.append((p[1].decode(), int(p[2])))
            except ValueError:
                pass
    return out

def coarse_index(url, total, step):
    offs = list(range(0, total, step))
    res = [None] * len(offs)
    def job(i):
        o = offs[i]
        for _ in range(4):
            try:
                pr = probe(url, o)
                res[i] = pr[0][0] if pr else '?'
                return
            except Exception:
                time.sleep(1)
        res[i] = '!'
    with ThreadPoolExecutor(8) as ex:
        list(ex.map(job, range(len(offs))))
    return list(zip(offs, res))

def first_off_with_chr(url, lo, hi, target):
    """在 [lo,hi] 字节区间二分，求首次出现 chr==target 的偏移"""
    while lo < hi:
        mid = (lo + hi) // 2
        pr = None
        for _ in range(5):
            try:
                pr = probe(url, mid); break
            except Exception:
                time.sleep(1)
        if pr is None:
            raise RuntimeError('probe fail %d' % mid)
        if pr and pr[0][0] == target:
            hi = mid
        else:
            lo = mid + 1
    return lo

def first_off_with_bp(url, lo, hi, target, bp):
    """在 [lo,hi] 区间二分，求该染色体上首个 BP >= bp 的偏移"""
    while lo < hi:
        mid = (lo + hi) // 2
        pr = None
        for _ in range(5):
            try:
                pr = probe(url, mid); break
            except Exception:
                time.sleep(1)
        if pr is None:
            raise RuntimeError('probe fail %d' % mid)
        ok = pr and pr[0][0] == target and pr[0][1] >= bp
        if ok:
            hi = mid
        else:
            lo = mid + 1
    return lo

if __name__ == '__main__':
    URL = sys.argv[1] if len(sys.argv) > 1 else \
        'https://www.columbiamedicine.org/divisions/kiryluk/gwas/IgA/META.IGA.LEVELS.ALL.COMBINED.txt'
    TOTAL = int(sys.argv[2]) if len(sys.argv) > 2 else 596831842
    STEP = 4 * 1024 * 1024
    t0 = time.time()
    idx = coarse_index(URL, TOTAL, STEP)
    print('coarse index: %d probes in %.1fs' % (len(idx), time.time() - t0))
    # 打印染色体边界
    prev = None
    for off, c in idx:
        if c != prev:
            print('  @%12d  %s' % (off, c))
            prev = c
    json.dump(idx, open('data/_chr_index.json', 'w'))
