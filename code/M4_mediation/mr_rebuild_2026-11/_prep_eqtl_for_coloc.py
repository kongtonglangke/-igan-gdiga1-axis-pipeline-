# -*- coding: utf-8 -*-
"""
_prep_eqtl_for_coloc.py
-----------------------
Normalise the eQTL Catalogue region cache into a plain 8-column table so that the
R colocalisation scripts can read it with fread().

This script does no statistics and no merging: it keeps only the five axis genes
and rewrites columns by position. Everything downstream (allele alignment,
coloc.abf) happens in R.

The cache files are long-format with an irregular first line and a truncated last
line, so rows are kept only when they carry the full 19 fields.

Write: data/eqtl_for_coloc/<QTD>.<chr>.tsv
       rsid  pos  ref  alt  maf  beta  se  gene
"""
import io
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

RAW = r"F:/longxiawenjian/IgAshenbing/00_rawdata/eQTL_Catalogue/region_cache"
OUT = r"F:/longxiawenjian/IgAshenbing/_mr_rebuild_20261110/data/eqtl_for_coloc"

ENSG = {"ENSG00000106392", "ENSG00000171155", "ENSG00000119514",
        "ENSG00000143641", "ENSG00000070731"}
QTDS = ["QTD000356", "QTD000261", "QTD000608", "QTD000607", "QTD000606"]
CHRS = ["1", "7", "9", "17", "X"]


def main():
    os.makedirs(OUT, exist_ok=True)
    total = 0
    for q in QTDS:
        for c in CHRS:
            src = os.path.join(RAW, "%s.%s.tsv" % (q, c))
            if not os.path.exists(src):
                continue
            dst = os.path.join(OUT, "%s.%s.tsv" % (q, c))
            kept = 0
            with io.open(src, encoding="utf-8", errors="replace") as fi, \
                    io.open(dst, "w", encoding="utf-8", newline="\n") as fo:
                fo.write("rsid\tpos\tref\talt\tmaf\tbeta\tse\tgene\n")
                for line in fi:
                    fld = line.rstrip("\n").split("\t")
                    if len(fld) != 19:
                        continue
                    if fld[0] not in ENSG:
                        continue
                    fo.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                        fld[18], fld[2], fld[3], fld[4], fld[7], fld[9], fld[10], fld[0]))
                    kept += 1
            total += kept
            print("  %s.%s  %8d rows" % (q, c, kept))
    print("total rows kept:", total)


if __name__ == "__main__":
    main()
