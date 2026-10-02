# -*- coding: utf-8 -*-
"""从本地已下载分块（data/META_IGA_LEVELS.txt.blk00000..）提取 chr1(GALNT2)/chr17(ST6GALNAC2) 区段。"""
import os, glob, sys

D = 'data'
OUT = 'data/igaLevels_local_part.tsv'
GENES = {
    'GALNT2':     ('1',  230193536, 230417870),
    'ST6GALNAC2': ('17', 74559792, 74583038),
}
PAD = 500000
want = {}
for g, (c, s, e) in GENES.items():
    want.setdefault(c, []).append((s - PAD, e + PAD, g))

blks = sorted(glob.glob(os.path.join(D, 'META_IGA_LEVELS.txt.blk*')))
print('local blocks:', len(blks))
rows = {g: [] for g in GENES}
for fp in blks:
    with open(fp, 'rb') as f:
        for ln in f:
            p = ln.rstrip(b'\n').split(b'\t')
            if len(p) < 8:
                continue
            try:
                c = p[1].decode().strip()
                bp = int(p[2])
            except ValueError:
                continue
            for (lo, hi, g) in want.get(c, []):
                if lo <= bp <= hi:
                    rows[g].append(ln if ln.endswith(b'\n') else ln + b'\n')
                    break
with open(OUT, 'wb') as g:
    g.write(b'SNP\tCHR\tBP_hg19\tA1\tA2\tBETA\tSE\tP\n')
    for k in ('GALNT2', 'ST6GALNAC2'):
        for r in rows[k]:
            g.write(r)
        print('  %-12s n=%d' % (k, len(rows[k])))
print('saved ->', OUT, os.path.getsize(OUT), 'bytes')
