# -*- coding: utf-8 -*-
# Patch X (2026-09-11): P1-2 — add axis_composite_null_q95 / axis_composite_emp_p to Table S10.
# Source: reproducibility/results/sc_regulatory/regulon_axis_composite.tsv (950 rows, target=AxisComposite).
# Values are filled ONLY on delivered rows with target=AxisComposite; per-gene rows stay empty
# (the empirical null belongs to the composite-axis test, not to per-gene tests).
# Also updates the A1 index wording for Table S10. Original S10 archived (留档).
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import csv, shutil, sys
from pathlib import Path

BASE = _paths.LEGACY
SRC = BASE("reproducibility/results/sc_regulatory/regulon_axis_composite.tsv")
DST = BASE("submission/additional_file_1/Table_S10_regulon_axis_associations.tsv")
ARC = BASE("submission/additional_file_1/Table_S10_regulon_axis_associations.before_nullq95_2026-09-11.tsv")
IDX = BASE("submission/additional_file_1/Additional_File_1_Figure_Legends.md")

# ---- 1. read source composite stats ----
with SRC.open(encoding="utf-8", newline="") as f:
    src = list(csv.DictReader(f, delimiter="\t"))
assert len(src) == 950, f"source rows {len(src)} != 950"
assert all(r["target"] == "AxisComposite" for r in src), "source has non-AxisComposite target"
smap = {}
for r in src:
    k = (r["b_subset"], r["tf"])
    assert k not in smap, f"dup source key {k}"
    smap[k] = (r["null_q95"], r["emp_p"])

# headline spot-checks from the manuscript (plasma DNMT1 0.025; naive KLF2 0.006; SP1 0.019)
def get(sub, tf):
    return smap[(sub, tf)]
assert get("B_plasma", "DNMT1")[1] == "0.025", get("B_plasma", "DNMT1")
assert get("B_naive", "KLF2")[1] == "0.006", get("B_naive", "KLF2")
assert get("B_naive", "SP1")[1] == "0.019", get("B_naive", "SP1")

# ---- 2. read delivered S10 ----
with DST.open(encoding="utf-8", newline="") as f:
    rdr = csv.DictReader(f, delimiter="\t")
    rows = list(rdr)
    cols = list(rdr.fieldnames)
assert len(rows) == 5700, f"S10 rows {len(rows)} != 5700"
assert cols == ["b_subset", "n_donors", "tf", "target", "rho", "p", "n_targets", "q"], cols
comp = [r for r in rows if r["target"] == "AxisComposite"]
assert len(comp) == 950, f"AxisComposite rows {len(comp)} != 950"
assert {(r["b_subset"], r["tf"]) for r in comp} == set(smap), "key mismatch S10 vs source"

# ---- 3. archive + merge ----
if not ARC.exists():
    shutil.copy2(DST, ARC)
newcols = cols + ["axis_composite_null_q95", "axis_composite_emp_p"]
n_fill = 0
for r in rows:
    if r["target"] == "AxisComposite":
        nq, ep = smap[(r["b_subset"], r["tf"])]
        r["axis_composite_null_q95"] = nq
        r["axis_composite_emp_p"] = ep
        n_fill += 1
    else:
        r["axis_composite_null_q95"] = ""
        r["axis_composite_emp_p"] = ""
assert n_fill == 950
with DST.open("w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=newcols, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)

# verify round-trip: DNMT1 plasma row carries 0.8810256... / 0.025
with DST.open(encoding="utf-8", newline="") as f:
    chk = [r for r in csv.DictReader(f, delimiter="\t")
           if r["b_subset"] == "B_plasma" and r["tf"] == "DNMT1" and r["target"] == "AxisComposite"]
assert len(chk) == 1 and chk[0]["axis_composite_emp_p"] == "0.025", chk
print(f"S10 merge OK: {n_fill} composite rows annotated; archive at {ARC.name}")

# ---- 4. A1 index wording ----
t = IDX.read_text(encoding="utf-8")
old = "| Table S10 | CollecTRI regulator-by-B-cell-state associations with each axis-gene target and the composite axis score (Spearman \u03c1, nominal P, Benjamini\u2013Hochberg q and number of panel targets; 950 regulator\u2013state combinations \u00d7 six targets = 5,700 rows). The target-matched random-set null and the empirical P values are summarised in the Fig. S6 legend and in the main text. Corresponds to Fig. 1c and Additional file 1: Fig. S6. |"
new = "| Table S10 | CollecTRI regulator-by-B-cell-state associations with each axis-gene target and the composite axis score (Spearman \u03c1, nominal P, Benjamini\u2013Hochberg q and number of panel targets; 950 regulator\u2013state combinations \u00d7 six targets = 5,700 rows). Composite-axis rows additionally carry the target-matched random-set null 95th percentile (axis_composite_null_q95) and the empirical P value (axis_composite_emp_p) from 1,000 matched random gene sets; per-gene rows leave these two fields empty because the empirical null applies to the composite-axis test. Corresponds to Fig. 1c and Additional file 1: Fig. S6. |"
assert t.count(old) == 1, "index S10 row not found verbatim"
IDX.write_text(t.replace(old, new), encoding="utf-8")
print("A1 index S10 wording OK")
print("ALL PATCH-X EDITS OK")
