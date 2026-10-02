# -*- coding: utf-8 -*-
"""从 WBBC VCF.gz 提取目标位点的中国人群等位频率（含 North/Central/South 亚群）。
用法: python extract_wbbc.py <vcf.gz> <chr:pos_hg38> [<chr:pos_hg38> ...]
输出 TSV: chrom  pos  ref  alt  AF  AN  NS  North_AF  Central_AF  South_AF  (原始INFO)
"""
import sys, gzip

def main():
    path = sys.argv[1]
    targets = {}
    for t in sys.argv[2:]:
        c, p = t.split(':')
        targets.setdefault(c, set()).add(int(p))
    print('targets:', targets, flush=True)
    hits = []
    n = 0
    with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('#'):
                continue
            n += 1
            c = line[:line.find('\t')]
            if c.lower().startswith('chr'):
                c = c[3:]
            if c in targets:
                i = line.find('\t'); j = line.find('\t', i + 1)
                try:
                    pos = int(line[i + 1:j])
                except ValueError:
                    continue
                if pos in targets[c]:
                    p = line.rstrip('\n').split('\t')
                    info = p[7] if len(p) > 7 else ''
                    d = {}
                    for kv in info.split(';'):
                        if '=' in kv:
                            k, v = kv.split('=', 1); d[k] = v
                    hits.append([c, str(pos), p[3], p[4], d.get('AF', ''), d.get('AN', ''), d.get('NS', ''),
                                 d.get('North_AF', ''), d.get('Central_AF', ''), d.get('South_AF', '')])
                    print('HIT', hits[-1], flush=True)
    print('scanned_lines=%d hits=%d' % (n, len(hits)), flush=True)
    with open(path + '.hits.tsv', 'w', encoding='utf-8') as g:
        g.write('CHROM\tPOS_hg38\tREF\tALT\tAF\tAN\tNS\tNorth_AF\tCentral_AF\tSouth_AF\n')
        for h in hits:
            g.write('\t'.join(h) + '\n')

if __name__ == '__main__':
    main()
