# -*- coding: utf-8 -*-
"""Recon: dump headers/shapes of all underlying data files for audit r2."""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import os, json, glob

ROOT = _paths.LEGACY
SUB = ROOT(r"阶段4.5_复现包\submission\additional_file_1")
DATA = ROOT(r"阶段4.5_复现包\reproducibility\data")
RES = ROOT(r"阶段4.5_复现包\reproducibility\results")
M1 = ROOT(r"阶段1\M1_GEO表达锚点")

def head(path, n=3):
    print("="*80)
    print(path)
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            lines = f.read().splitlines()
        print(f"[{len(lines)} lines]")
        for ln in lines[:n]:
            print(ln[:600])
    except Exception as e:
        print("ERR", e)

# submission supp tables
for f in ["Table_S1_eQTL_5x5_matrix.tsv","Table_S2a_coloc_25comparisons.tsv",
          "Table_S2b_coloc_verdict_14computable.tsv",
          "Table_S5a_sigCpG_5CpG_pos.tsv","Table_S5b_MR_CpG_x_IgAN.tsv",
          "Table_S5c_coloc_mQTL_x_IgAN.tsv","Table_S5d_coloc_mQTL_x_eQTL.tsv",
          "Table_S6a_intervenability_scores.tsv","Table_S6b_drug_directory.tsv",
          "Table_S11_Bcell_eQTL_heterogeneity.tsv","Table_S12_EAS_eQTL_crossref.tsv",
          "Table_S13_lead_C1GALT1_eQTL_EAS.tsv"]:
    head(os.path.join(SUB, f), 4)

# data dir tree
print("#"*80)
for dp, dn, fn in os.walk(DATA):
    for f in sorted(fn):
        print(os.path.join(dp, f).replace(DATA, "DATA"))

# M1 DE source headers
for f in ["GSE73953_IgAN_vs_HC_allprobe.csv","GSE115857_IgAN_vs_Ctrl_LD_allprobe.csv",
          "GSE93798_IgAN_vs_Ctrl_allprobe.csv"]:
    head(os.path.join(M1, f), 2)

# key JSONs
for f in [r"sc_regulatory\step1_summary.json", r"sc_regulatory\step2_summary.json",
          r"sc_regulatory\step3_summary.json", r"B2_bcell_eqtl\step1_summary.json",
          r"C1_eas_eqtl\step1_summary.json", r"C1_eas_eqtl\step2_summary.json"]:
    p = os.path.join(RES, f)
    print("="*80); print(p)
    try:
        with open(p, encoding="utf-8") as fh:
            print(fh.read()[:2500])
    except Exception as e:
        print("ERR", e)
