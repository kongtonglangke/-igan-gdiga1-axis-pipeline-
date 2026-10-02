# -*- coding: utf-8 -*-
"""小块动态分片下载器 v2：失败重试只重下小块，支持断点续传（已完成的 .blk 跳过）。
用法: python dl_par2.py <url> <out_path> [nthreads] [blocksize]
"""
import sys, os, ssl, time, queue, threading, urllib.request
from concurrent.futures import ThreadPoolExecutor

CTX = ssl.create_default_context(); CTX.check_hostname = False; CTX.verify_mode = ssl.CERT_NONE
UA = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'}

def log(m):
    print(time.strftime('[%H:%M:%S] ') + str(m), flush=True)

def size_of(url):
    r = urllib.request.Request(url, method='HEAD', headers=UA)
    with urllib.request.urlopen(r, timeout=90, context=CTX) as f:
        return int(f.headers['Content-Length'])

def fetch(url, a, b, fp, tries=12):
    need = b - a + 1
    for t in range(tries):
        try:
            h = dict(UA); h['Range'] = 'bytes=%d-%d' % (a, b)
            r = urllib.request.Request(url, headers=h)
            with urllib.request.urlopen(r, timeout=120, context=CTX) as f, open(fp, 'wb') as g:
                while True:
                    c = f.read(131072)
                    if not c:
                        break
                    g.write(c)
            if os.path.getsize(fp) == need:
                return True
        except Exception:
            pass
        time.sleep(1.5 + t * 0.5)
    return False

def main():
    url, out = sys.argv[1], sys.argv[2]
    nthreads = int(sys.argv[3]) if len(sys.argv) > 3 else 16
    bs = int(sys.argv[4]) if len(sys.argv) > 4 else 1_000_000
    total = size_of(url)
    log('size=%d (%.1f MB) threads=%d block=%d' % (total, total / 1048576.0, nthreads, bs))
    if os.path.exists(out) and os.path.getsize(out) == total:
        log('DONE already'); return
    jobs = []
    a, i = 0, 0
    while a < total:
        b = min(a + bs - 1, total - 1)
        jobs.append((i, a, b, out + '.blk%05d' % i)); a = b + 1; i += 1
    todo = [j for j in jobs if not (os.path.exists(j[3]) and os.path.getsize(j[3]) == j[2] - j[1] + 1)]
    log('blocks=%d todo=%d' % (len(jobs), len(todo)))
    q = queue.Queue()
    for j in todo:
        q.put(j)
    st = {'n': 0, 'fail': []}; lock = threading.Lock()
    def worker():
        while True:
            try:
                j = q.get_nowait()
            except queue.Empty:
                return
            ok = fetch(url, j[1], j[2], j[3])
            with lock:
                st['n'] += 1
                if not ok: st['fail'].append(j[0])
                if st['n'] % 40 == 0: log('  %d/%d' % (st['n'], len(todo)))
            if not ok: log('  FAIL blk %d' % j[0])
    with ThreadPoolExecutor(nthreads) as ex:
        list(ex.map(lambda _: worker(), range(nthreads)))
    if st['fail']:
        log('FAILED blocks: %s' % st['fail'][:20]); sys.exit(1)
    log('merging ...')
    with open(out, 'wb') as g:
        for _, _, _, fp in jobs:
            with open(fp, 'rb') as f:
                while True:
                    c = f.read(1 << 22)
                    if not c: break
                    g.write(c)
    if os.path.getsize(out) != total:
        log('SIZE MISMATCH'); sys.exit(1)
    for _, _, _, fp in jobs: os.remove(fp)
    log('DONE %s %d bytes' % (out, os.path.getsize(out)))

if __name__ == '__main__':
    main()
