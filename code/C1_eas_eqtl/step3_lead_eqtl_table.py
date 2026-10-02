#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
step3_lead_eqtl_table.py — C1 第三步：Gd-IgA1 lead 变异在 ImmuNexUT（东亚）中的 C1GALT1 eQTL 明细表
================================================================================
输入：results/C1_eas_eqtl/axis_hits_raw.tsv
输出：results/C1_eas_eqtl/TableS13_lead_C1GALT1_eQTL_EAS.tsv
      每 (lead × 细胞亚型) 一行，仅保留该 lead 作为 C1GALT1 条件独立信号出现的记录；
      另行为 "absent" 占位（若该 lead 在 B 系亚型中不是独立信号）。
"""

# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import os
from collections import defaultdict

BASE = _paths.LEGACY
RES = BASE("阶段4.5_复现包", "reproducibility", "results", "C1_eas_eqtl")

LEADS = {"rs10238682": ("chr7", 7215386), "rs13226913": ("chr7", 7207215)}
B_SUBSETS = ["Naive_B", "USM_B", "SM_B", "DN_B", "Plasmablast"]

recs = defaultdict(list)
all_subsets = set()

with open(os.path.join(RES, "axis_hits_raw.tsv"), encoding="utf-8") as fh:
    next(fh)
    for line in fh:
        member, _, payload = line.partition("\t")
        p = payload.rstrip("\n").split("\t")
        if len(p) < 14 or p[1] != "C1GALT1":
            continue
        sub = member.split("_conditional")[0]
        all_subsets.add(sub)
        for rs in LEADS:
            if rs in p[5]:
                recs[rs].append({
                    "lead": rs, "cell_subset": sub,
                    "rank_of_association": int(float(p[9])),
                    "nominal_P": float(p[10]),
                    "slope": float(p[11]),
                    "variant_id": p[5], "variant_pos": int(p[7]),
                })

rows = []
for rs in LEADS:
    hits = {r["cell_subset"]: r for r in recs[rs]}
    for sub in sorted(all_subsets, key=lambda s: (s not in B_SUBSETS, s)):
        if sub in hits:
            r = hits[sub]
            rows.append({
                "lead_rsid": rs, "cell_subset": sub, "B_lineage": "yes" if sub in B_SUBSETS else "no",
                "is_independent_signal": "yes", "rank_of_association": r["rank_of_association"],
                "nominal_P": r["nominal_P"], "effect_slope": r["slope"],
                "variant_pos_hg38": r["variant_pos"],
            })
        elif sub in B_SUBSETS:
            rows.append({
                "lead_rsid": rs, "cell_subset": sub, "B_lineage": "yes",
                "is_independent_signal": "no", "rank_of_association": "",
                "nominal_P": "", "effect_slope": "", "variant_pos_hg38": LEADS[rs][1],
            })

cols = ["lead_rsid", "cell_subset", "B_lineage", "is_independent_signal",
        "rank_of_association", "nominal_P", "effect_slope", "variant_pos_hg38"]
out = os.path.join(RES, "TableS13_lead_C1GALT1_eQTL_EAS.tsv")
with open(out, "w", encoding="utf-8") as f:
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join(str(r[c]) for c in cols) + "\n")

print(f"[done] {out}  rows={len(rows)}")
print("\nB-lineage summary:")
for rs in LEADS:
    b = [r for r in rows if r["lead_rsid"] == rs and r["B_lineage"] == "yes"]
    hit = [r for r in b if r["is_independent_signal"] == "yes"]
    print(f"  {rs}: {len(hit)}/{len(b)} B-lineage subsets are independent signals")
    for r in hit:
        print(f"     {r['cell_subset']:12s} rank={r['rank_of_association']} P={r['nominal_P']:.3g} slope={r['effect_slope']:+.3f}")
