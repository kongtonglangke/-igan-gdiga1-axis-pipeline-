# -*- coding: utf-8 -*-
"""轴基因区增量分析：C 项（亚洲 IgAN GWAS 位点直查）+ D 项（血清 IgA GWAS 轴区扫描）。
用法: python analyze_axis.py <IgAN_Asian_only.txt> <igaLevels_axisChrs.tsv> <out_dir>
"""
import sys, os, json

# 4 个已发表 Gd-IgA1 lead（hg19，与两文件一致）
LEADS = [
    ('rs13226913', '7', 7246846),
    ('rs10238682', '7', 7255017),
    ('rs5910940',  'X', 119758810),
    ('rs7856182',  '9', 101633312),
]
# 五基因区 hg19（GRCh37，含 ±500 kb 上下文）
GENES = {
    'C1GALT1':    ('7',  7196565, 7288282),
    'C1GALT1C1':  ('X',  119759648, 119764005),
    'GALNT12':    ('9',  101569981, 101612363),
    'GALNT2':     ('1',  230193536, 230417870),
    'ST6GALNAC2': ('17', 74559792, 74583038),
}
PAD = 500000

def norm(c):
    c = str(c).strip()
    if c.lower().startswith('chr'):
        c = c[3:]
    if c == '23': c = 'X'
    return c

def scan_gwas(path, want):
    """want: dict chr -> list of (lo,hi,label)。返回 {label: [(chr,bp,a1,a2,beta,se,p)]}"""
    res = {lbl: [] for lst in want.values() for (_, _, lbl) in lst}
    if not os.path.exists(path):
        print('  !! missing', path); return res
    with open(path, encoding='utf-8', errors='replace') as f:
        head = f.readline()
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) < 8: continue
            c = norm(p[1])
            try: bp = int(p[2])
            except ValueError: continue
            for (lo, hi, lbl) in want.get(c, []):
                if lo <= bp <= hi:
                    res[lbl].append((c, bp, p[3], p[4], p[5], p[6], p[7]))
    return res

def main():
    asian, igalev, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
    os.makedirs(outdir, exist_ok=True)
    report = {}
    # ---------- C 项 ----------
    print('=' * 70); print('C 项：亚洲 IgAN GWAS (Asian-only meta) 直查 4 个 lead')
    want = {}
    for rs, c, bp in LEADS:
        want.setdefault(c, []).append((bp, bp, rs))
    cres = scan_gwas(asian, want)
    ctab = []
    for rs, c, bp in LEADS:
        rows = cres.get(rs, [])
        if rows:
            cc, bb, a1, a2, beta, se, p = rows[0]
            ctab.append((rs, cc, bb, a1, a2, beta, se, p))
            print('  %-12s %s:%d  A1=%s A2=%s  BETA=%s SE=%s  P=%s' % (rs, cc, bb, a1, a2, beta, se, p))
        else:
            ctab.append((rs, c, bp, '-', '-', '-', '-', 'absent'))
            print('  %-12s %s:%d  -> 不在文件中（absent）' % (rs, c, bp))
    report['C_leads'] = ctab
    # ---------- D 项 ----------
    print('=' * 70); print('D 项：血清总 IgA GWAS (Liu 2022, n=41,263) 轴基因区 ±500kb 扫描')
    want2 = {}
    for g, (c, s, e) in GENES.items():
        want2.setdefault(c, []).append((s - PAD, e + PAD, g))
    dres = scan_gwas(igalev, want2)
    dtab = []
    for g, (c, s, e) in GENES.items():
        rows = dres.get(g, [])
        if rows:
            rows_sorted = sorted(rows, key=lambda r: float(r[6]) if r[6] not in ('NA', '') else 1)
            best = rows_sorted[0]
            dtab.append((g, c, s, e, len(rows), best[1], best[2], best[3], best[4], best[5], best[6]))
            print('  %-12s chr%s:%-10d-%-10d  nSNP=%-6d  minP @ %s:%d  A1=%s A2=%s BETA=%s SE=%s P=%s' %
                  (g, c, s, e, len(rows), best[0], best[1], best[2], best[3], best[4], best[5], best[6]))
        else:
            dtab.append((g, c, s, e, 0, '-', '-', '-', '-', '-', '-'))
            print('  %-12s chr%s:%-10d-%-10d  nSNP=0' % (g, c, s, e))
    report['D_axis_regions'] = dtab
    with open(os.path.join(outdir, 'analysis_axis_summary.json'), 'w', encoding='utf-8') as g:
        json.dump(report, g, indent=2, ensure_ascii=False)
    print('=' * 70); print('saved ->', os.path.join(outdir, 'analysis_axis_summary.json'))

if __name__ == '__main__':
    main()
