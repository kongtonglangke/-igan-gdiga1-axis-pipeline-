# -*- coding: utf-8 -*-
"""M7 step 5 - build the supplementary tables from the raw inputs.

Rebuilds `Table_S14_ancestry_matched_lookups.tsv` deterministically from:
  * the Asian-only IgAN summary statistics (Kiryluk lab; hg19; SNP key = CHR:BP)
  * the WBBC Chinese population allele-frequency hits  (chr7 / chr9 TSVs)
  * the gnomAD v4 EAS/NFE + Han Gd-IgA1 + IgAN-meta columns already archived
    under reproducibility/data/M5_eastasia/

Usage (from the analysis workspace root):
  python step_m7_5_build_tables.py [--asian-only PATH] [--m5 PATH]
                                   [--wbbc7 PATH] [--wbbc9 PATH] [--out PATH]

Formatting rules reproduce the delivered table exactly:
  gnomAD EAS/NFE                     %.3f
  WBBC Chinese (incl. sub-populations)   %.5f when AF < 0.1 else %.4f
  beta / SE                          verbatim from the summary-statistics file
  Asian-only P                       %.3g
  Han Gd-IgA1 P                      %.2g  (exponent normalised, e-09 -> e-9)
  European IgAN P                    %.3f  (NA when the variant is untestable)
"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import argparse
import os

LEADS = [
    # rsid,        gene,        chr, hg19 bp,   allele (gnomAD ALT)
    ('rs13226913', 'C1GALT1',   '7', 7246846,  'C'),
    ('rs10238682', 'C1GALT1',   '7', 7255017,  'G'),
    ('rs7856182',  'GALNT12',   '9', 101633312, 'T'),
    ('rs5910940',  'C1GALT1C1', 'X', 119758810, 'A'),
]
HEADER = ('variant\tgene\tallele\tgnomAD_v4_EAS_AF\tgnomAD_v4_NFE_AF\tWBBC_Chinese_AF\t'
          'WBBC_north_AF\tWBBC_central_AF\tWBBC_south_AF\tAsianOnlyIgAN_A1\tAsianOnlyIgAN_A2\t'
          'AsianOnlyIgAN_BETA\tAsianOnlyIgAN_SE\tAsianOnlyIgAN_P\tHanGdIgA1_P\tEuropeanIgAN_P')


def af(x):
    return ('%.5f' % x) if abs(x) < 0.1 else ('%.4f' % x)


def g2(x):
    s = '%.2g' % x
    return s.replace('e-0', 'e-').replace('e+0', 'e+')


def load_asian(path):
    want = {('%s:%d' % (c, bp)): rs for rs, _g, c, bp, _a in LEADS}
    out = {}
    with open(path, encoding='utf-8', errors='replace') as f:
        for i, line in enumerate(f):
            if i == 0 and line[:3] == 'SNP':
                continue
            q = line.rstrip('\n').split('\t')
            if len(q) < 8:
                continue
            if q[0] in want:
                out[want[q[0]]] = q            # SNP CHR BP A1 A2 BETA SE P
    return out


def load_gnomad(path):
    out = {}
    with open(path, encoding='utf-8') as f:
        hdr = f.readline().rstrip('\n').split('\t')
        idx = {k: i for i, k in enumerate(hdr)}
        for line in f:
            q = line.rstrip('\n').split('\t')
            if len(q) <= idx['IgAN_meta_P']:
                q += [''] * (idx['IgAN_meta_P'] + 1 - len(q))
            out[q[idx['lead_rsid']]] = q
    return out, idx


def load_wbbc(files):
    """rsid-agnostic: WBBC hits are keyed by hg38 position; match on position."""
    pos = {}
    for p in files:
        if not os.path.exists(p):
            continue
        with open(p, encoding='utf-8') as f:
            for line in f:
                q = line.rstrip('\n').split('\t')
                if not q[0].isdigit():
                    continue
                # CHROM POS_hg38 REF ALT AF AN NS North_AF Central_AF South_AF
                pos[(q[0], int(q[1]))] = q
    return pos


# hg38 positions of the four leads (needed to join with the WBBC hits)
HG38 = {'rs13226913': ('7', 7207215), 'rs10238682': ('7', 7215386),
        'rs7856182': ('9', 98871030), 'rs5910940': ('X', 120624955)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--asian-only', default='data/IgAN_Asian_only.txt')
    ap.add_argument('--m5', default=_paths.data(
        'M5_eastasia', 'M5_gnomAD_EASEUR_频率_20260903.tsv'))
    ap.add_argument('--wbbc7', default='data/WBBC.chr7.vcf.gz.hits.tsv')
    ap.add_argument('--wbbc9', default='data/WBBC.chr9.vcf.gz.hits.tsv')
    ap.add_argument('--out', default='out/Table_S14_ancestry_matched_lookups.tsv')
    a = ap.parse_args()

    asian = load_asian(a.asian_only)
    gno, gi = load_gnomad(a.m5)
    wbbc = load_wbbc([a.wbbc7, a.wbbc9])

    lines = [HEADER]
    for rs, gene, c, bp, allele in LEADS:
        g = gno[rs]
        eas = float(g[gi['gnomAD_EAS_AF']])
        nfe = float(g[gi['gnomAD_EUR_NFE_AF']])
        hanp = float(g[gi['Wang2021_Han_P']]) if g[gi['Wang2021_Han_P']] else None
        eurp = g[gi['IgAN_meta_P']]
        w = wbbc.get(HG38[rs])
        if w:
            wch, wn, wc, ws = (float(w[4]), float(w[7]), float(w[8]), float(w[9]))
            wf = [af(wch), af(wn), af(wc), af(ws)]
        else:
            wf = ['NA'] * 4
        z = asian.get(rs)
        if z:
            a1, a2, beta, se, p = z[3], z[4], z[5], z[6], float(z[7])
            ap_ = '%.3g' % p
        else:
            a1 = a2 = beta = se = ap_ = 'NA'
        lines.append('\t'.join([
            rs, gene, allele,
            '%.3f' % eas, '%.3f' % nfe,
            wf[0], wf[1], wf[2], wf[3],
            a1, a2, beta, se, ap_,
            (g2(hanp) if hanp is not None else 'NA'),
            ('%.3f' % float(eurp)) if eurp else 'NA',
        ]))

    data = ('\n'.join(lines) + '\n').encode('utf-8')
    outdir = os.path.dirname(a.out)
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    open(a.out, 'wb').write(data)
    print('wrote', a.out, len(data), 'B')
    print(data.decode('utf-8'))


if __name__ == '__main__':
    main()
