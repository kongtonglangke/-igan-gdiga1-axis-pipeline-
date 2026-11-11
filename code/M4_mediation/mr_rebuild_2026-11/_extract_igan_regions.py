# -*- coding: utf-8 -*-
"""
_extract_igan_regions.py
------------------------
One-off extraction of the four IgAN GWAS regions needed by the R colocalisation
scripts (C1GALT1 on chr7, GALNT12 on chr9, GALNT2 on chr1, ST6GALNAC2 on chr17;
C1GALT1C1 sits on chrX and GCST90018866 has no X chromosome, so it is skipped).

This step only moves bytes - all statistics are computed in R with coloc::coloc.abf.

Write: data/igan_regions_for_coloc/<chr>.tsv
"""
import gzip
import io
import os
import sys
import time

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

RAW = r"F:/longxiawenjian/IgAshenbing/00_rawdata/GCST90018866/GCST90018866.h.tsv.gz"
OUT = r"F:/longxiawenjian/IgAshenbing/_mr_rebuild_20261110/data/igan_regions_for_coloc"

# gene body +- 1 Mb, matching the eQTL side flank used in step1_coloc_abf.py
REGIONS = {
    "7":  (6_156_934, 8_259_333),      # C1GALT1
    "9":  (97_807_613, 99_851_267),    # GALNT12
    "1":  (229_057_990, 231_282_128),  # GALNT2
    "17": (75_564_424, 77_586_956),    # ST6GALNAC2
}


def main():
    os.makedirs(OUT, exist_ok=True)
    fh = {c: io.open(os.path.join(OUT, "chr%s.tsv" % c), "w", encoding="utf-8", newline="\n")
          for c in REGIONS}
    for f in fh.values():
        f.write("chr\tpos\tea\toa\tbeta\tse\teaf\tp\tvid\thm_conv\thm_code\trsid\n")

    n = kept = 0
    t0 = time.time()
    with gzip.open(RAW, "rt", encoding="utf-8", errors="replace") as f:
        f.readline()
        for line in f:
            n += 1
            fld = line.rstrip("\n").split("\t")
            if len(fld) < 12:
                continue
            c = fld[0]
            if c not in REGIONS:
                continue
            try:
                pos = int(fld[1])
            except ValueError:
                continue
            lo, hi = REGIONS[c]
            if lo <= pos <= hi:
                fh[c].write("\t".join(fld) + "\n")
                kept += 1
            if n % 5_000_000 == 0:
                print("  scanned %.0fM lines, kept %d, %.0fs" % (n / 1e6, kept, time.time() - t0),
                      flush=True)

    for f in fh.values():
        f.close()
    print("scanned %d lines, kept %d variants, %.0fs" % (n, kept, time.time() - t0))
    for c in sorted(REGIONS):
        p = os.path.join(OUT, "chr%s.tsv" % c)
        with io.open(p, encoding="utf-8") as f:
            print("  chr%-2s %8d rows  %s" % (c, sum(1 for _ in f) - 1, p))


if __name__ == "__main__":
    main()
