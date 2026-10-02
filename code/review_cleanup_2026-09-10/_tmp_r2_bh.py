# -*- coding: utf-8 -*-
"""BH sensitivity recheck over the 13-lookup family (S4a leads + S4b windows)."""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import csv, os, sys

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_tmp_r2_bh_out.txt")
sys.stdout = open(OUT, "w", encoding="utf-8")

BASE = _paths.at("阶段4.5_复现包/submission/additional_file_1")
P = []
# 3 evaluable leads from M2b lead_x_IgAN_Pvalue.tsv
LEAD = _paths.legacy("阶段4.5_复现包/reproducibility/data/M2b_susie_coloc/lead_x_IgAN_Pvalue.tsv")
with open(LEAD, encoding="utf-8") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["igan_p"]:
            P.append((float(r["igan_p"]), "lead:" + r["lead_rsid"]))
with open(BASE("Table_S4b_cisQTL_topvar_IgAN_locus_lookup.tsv"), encoding="utf-8") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        P.append((float(r["locus_minP_P"]), "S4b:" + r["gene"]))

P.sort(key=lambda x: x[0])
m = len(P)
print("m =", m)
raw = [p * m / i for i, (p, _) in enumerate(P, 1)]
q = [min(raw[i - 1:]) for i in range(1, m + 1)]
for i, ((p, tag), qq) in enumerate(zip(P, q), 1):
    print(f"rank{i:2d} {tag:28s} P={p:.4e}  Px13/{i}={raw[i-1]:.4e}  q_stepup={qq:.4e}")
print("=> smallest window q:", min(q))
