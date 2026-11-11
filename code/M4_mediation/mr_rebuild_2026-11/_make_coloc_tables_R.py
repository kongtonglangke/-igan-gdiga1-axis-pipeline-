# -*- coding: utf-8 -*-
"""
_make_coloc_tables_R.py
-----------------------
Turn the two R colocalisation outputs (coloc::coloc.abf) into the exact file
formats the manuscript, the figures and the supplement already use, so that every
region-level colocalisation in the package comes from the official coloc package.

  13_coloc_eqtl_igan_R.tsv  ->  coloc_5genes_x_contexts_R.tsv      (25 rows, Fig. 3B/C, S2a)
                            ->  Table_S19_..._R.tsv                (14 computable, 3 dp)
  14_coloc_mqtl_eqtl_R.tsv  ->  Table_S5d_coloc_mQTL_x_eQTL_R.tsv  (25 rows)

Old Python outputs are left untouched; the new files carry the _R suffix.
"""
import csv
import io
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

MR = r"F:/longxiawenjian/IgAshenbing/_mr_rebuild_20261110"
STAGE2 = r"F:/longxiawenjian/IgAshenbing/阶段2/主证据2_共享结构"
OUT = os.path.join(MR, "out", "supp_tables")


def read(p):
    with io.open(p, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main():
    os.makedirs(OUT, exist_ok=True)

    # ---------- eQTL x IgAN ----------
    a = read(os.path.join(MR, "out", "13_coloc_eqtl_igan_R.tsv"))
    a.sort(key=lambda r: (r["gene"], r["qtd"]))
    cols = ["gene", "qtd", "context", "n_SNP", "PPH0", "PPH1", "PPH2", "PPH3", "PPH4", "note"]
    p1 = os.path.join(STAGE2, "coloc_5genes_x_contexts_R.tsv")
    with io.open(p1, "w", encoding="utf-8", newline="\n") as f:
        f.write("\t".join(cols) + "\n")
        for r in a:
            f.write("\t".join([r["gene"], r["qtd"], r["context"], r["n_snp"],
                               r["PPH0"], r["PPH1"], r["PPH2"], r["PPH3"], r["PPH4"],
                               r["note"]]) + "\n")

    # ---------- S19, 14 computable, 3 dp ----------
    p2 = os.path.join(OUT, "Table_S19_coloc_eQTL_x_IgAN_R.tsv")
    ctx = {"GTEx_v10_blood": "GTEx v10 whole blood",
           "GTEx_v10_kidney_cortex": "GTEx v10 kidney cortex",
           "OneK1K_B_naive": "OneK1K naive B",
           "OneK1K_B_memory": "OneK1K memory B",
           "OneK1K_B_intermediate": "OneK1K intermediate B"}
    kept = [r for r in a if r["PPH4"] not in ("", "NA")]
    with io.open(p2, "w", encoding="utf-8", newline="\n") as f:
        f.write("Gene\tContext\tn_SNP\tPPH1\tPPH4\n")
        for r in kept:
            f.write("%s\t%s\t%s\t%.3f\t%.3f\n" % (
                r["gene"], ctx[r["context"]], r["n_snp"],
                float(r["PPH1"]), float(r["PPH4"])))

    # ---------- mQTL x eQTL ----------
    b = read(os.path.join(MR, "out", "14_coloc_mqtl_eqtl_R.tsv"))
    b.sort(key=lambda r: (r["cpg"], r["qtd"]))
    p3 = os.path.join(OUT, "Table_S5d_coloc_mQTL_x_eQTL_R.tsv")
    cols3 = ["cpg", "qtd", "context", "n_shared", "mqtl_top", "eqtl_top", "eqtl_top_P",
             "PPH0", "PPH1", "PPH2", "PPH3", "PPH4", "note"]
    with io.open(p3, "w", encoding="utf-8", newline="\n") as f:
        f.write("\t".join(cols3) + "\n")
        for r in b:
            f.write("\t".join([r[c] for c in cols3]) + "\n")

    print("25 eQTL x IgAN comparisons ->", p1)
    print("14 computable ->", p2)
    print("25 mQTL x eQTL comparisons ->", p3)
    print("\nR summary, eQTL x IgAN: max PPH4 = %.4f, PPH1 range %.3f-%.3f (14 computable)" % (
        max(float(r["PPH4"]) for r in kept),
        min(float(r["PPH1"]) for r in kept), max(float(r["PPH1"]) for r in kept)))
    print("R summary, mQTL x eQTL: %d comparisons >= 0.8" % sum(1 for r in b if float(r["PPH4"]) >= 0.8))
    for r in b:
        if float(r["PPH4"]) >= 0.8:
            print("   %-11s %-24s %.4f" % (r["cpg"], r["context"], float(r["PPH4"])))


if __name__ == "__main__":
    main()
