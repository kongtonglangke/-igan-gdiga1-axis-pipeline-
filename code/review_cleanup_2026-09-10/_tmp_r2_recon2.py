# -*- coding: utf-8 -*-
"""Recon2: full dump of small key files + grep manuscript numbers."""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import os, json

ROOT = _paths.LEGACY
SUB = ROOT(r"阶段4.5_复现包\submission\additional_file_1")
DATA = ROOT(r"阶段4.5_复现包\reproducibility\data")
RES = ROOT(r"阶段4.5_复现包\reproducibility\results")
M1 = ROOT(r"阶段1\M1_GEO表达锚点")
MS = ROOT(r"阶段4.5_复现包\submission\manuscript")

def dump(path, label=None, maxlines=40):
    print("="*90)
    print(label or path)
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            lines = f.read().splitlines()
        print(f"[{len(lines)} lines]")
        for ln in lines[:maxlines]:
            print(ln[:500])
    except Exception as e:
        print("ERR", e)

dump(os.path.join(SUB, "Table_S1_eQTL_5x5_matrix.tsv"), "S1 full")
dump(os.path.join(DATA, r"M2_lead_eqtl\lead_eQTL_matrix.tsv"), "lead_eQTL_matrix full")
dump(os.path.join(M1, "GSE73953_IgAN_vs_HC_best.csv"), "GSE73953 best")
dump(os.path.join(M1, "GSE115857_IgAN_vs_Ctrl_LD_best.csv"), "GSE115857 best")
dump(os.path.join(M1, "GSE93798_IgAN_vs_Ctrl_best.csv"), "GSE93798 best")
dump(os.path.join(DATA, r"M4_mediation\lead_x_mQTL_sigCpG_raw.tsv"), "M4 sigCpG raw")
dump(os.path.join(SUB, "Table_S6a_intervenability_scores.tsv"), "S6a full")
dump(os.path.join(DATA, r"M5_eastasia\M5_gnomAD_EASEUR_频率_20260903.tsv"), "M5 gnomAD")
dump(os.path.join(DATA, r"M5_eastasia\lead_x_IgAN_eaf.tsv"), "M5 IgAN eaf")
dump(os.path.join(DATA, r"M2b_susie_coloc\lead_x_Wang2021_GdIgA1_P.tsv"), "Wang2021 P")
dump(os.path.join(DATA, r"M2b_susie_coloc\lead_x_IgAN_Pvalue.tsv"), "IgAN meta P")
dump(os.path.join(RES, r"directional_validation\bcell_subset_counts.tsv"), "DV subset counts")

# step2_summary empirical P section
with open(os.path.join(RES, r"sc_regulatory\step2_summary.json"), encoding="utf-8") as f:
    js = json.load(f)
print("="*90)
print("step2_summary keys:", list(js.keys()))
for k in js:
    if k not in ("top_composite",):
        print(k, "->", json.dumps(js[k], ensure_ascii=False)[:800])
# find empirical P entries in top_composite or elsewhere
s = json.dumps(js, ensure_ascii=False)
for tf in ["DNMT1", "KLF2", "SP1", "MYC", "TFAP2A"]:
    idx = s.find(tf)
    while idx != -1:
        print(f"...{tf} ctx:", s[max(0,idx-120):idx+200])
        idx = s.find(tf, idx+1)
        break
