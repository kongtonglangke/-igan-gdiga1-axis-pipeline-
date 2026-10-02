#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
step2_summarise_crossref.py — C1 第二步：汇总 ImmuNexUT 轴基因 eQTL 交叉参照表
================================================================================
输入：results/C1_eas_eqtl/axis_hits_raw.tsv（step1 流式抽出的轴基因行）
      （列：member \t 制表符分隔的原始行）
      原始列（来自 E-GEAD-398 conditional_eQTL_FDR0.05.tar，QTLtools conditional pass）：
        1 Gene_id  2 Gene_name  3 CHR  4 TSS_position  5 Number_of_variants_cis
        6 Variant_ID  7 Variant_CHR  8 Variant_position_start  9 Variant_position_end
        10 Rank_of_association  11 Forward_nominal_P  12 Forward_slope
        13 Backward_P  14 Backward_slope
      注意：文件只含该细胞亚型中 FDR<0.05 的 eGene ⇒ 出现即该基因为该亚型的 eGene。

输出：TableS12_EAS_eQTL_crossref.tsv   每 (亚型 × 轴基因) 一行：独立信号数、最优变异、P、slope
      step2_summary.json              跨 28 亚型的 eGene 广度 + B 系亚型结论 + lead 位点查询
"""

# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import os
import json
from collections import defaultdict

BASE = _paths.LEGACY
RES = BASE("阶段4.5_复现包", "reproducibility", "results", "C1_eas_eqtl")

AXIS = {
    "ENSG00000106392": "C1GALT1",
    "ENSG00000171155": "C1GALT1C1",
    "ENSG00000143641": "GALNT2",
    "ENSG00000119514": "GALNT12",
    "ENSG00000070731": "ST6GALNAC2",
}
B_SUBSETS = ["Naive_B", "USM_B", "SM_B", "DN_B", "Plasmablast"]
# 已发表 Gd-IgA1 lead 变异的 hg38 位置（与 B2 脚本一致）
LEAD_POS = {"rs13226913": ("chr7", 7207215), "rs10238682": ("chr7", 7215386)}

rows = defaultdict(list)
subset_hits = defaultdict(int)

with open(os.path.join(RES, "axis_hits_raw.tsv"), encoding="utf-8") as fh:
    next(fh)  # header member\tline
    for line in fh:
        member, _, payload = line.partition("\t")
        p = payload.rstrip("\n").split("\t")
        if len(p) < 14:
            continue
        subset = member.split("_conditional")[0]
        ensg = p[0].split(".")[0]
        gene = p[1]
        if ensg not in AXIS:
            continue
        subset_hits[subset] += 1
        rows[(subset, AXIS[ensg])].append({
            "gene": AXIS[ensg],
            "ensg": ensg,
            "subset": subset,
            "tss": int(p[3]),
            "variant_id": p[5],
            "vchr": p[6],
            "vpos": int(p[7]),
            "rank": int(float(p[9])),
            "p": float(p[10]),
            "slope": float(p[11]),
            "bwd_p": float(p[12]),
        })

# ---- 每 (亚型×轴基因) 汇总 ----
out = []
for (subset, gene), recs in sorted(rows.items()):
    top = min(recs, key=lambda r: (r["rank"], r["p"]))
    best_p = min(r["p"] for r in recs)
    out.append({
        "cell_subset": subset,
        "gene": gene,
        "is_eGene_FDR0.05": "yes",
        "n_independent_signals": len(recs),
        "top_variant_rank0": top["variant_id"],
        "top_variant_pos": top["vpos"],
        "top_nominal_P": top["p"],
        "top_slope": top["slope"],
        "min_nominal_P_any_signal": best_p,
        "B_lineage": "yes" if subset in B_SUBSETS else "no",
    })

# ---- lead 位点是否出现在 C1GALT1 的独立信号中 ----
lead_query = []
c1 = [r for r in sum(rows.values(), []) if r["gene"] == "C1GALT1"]
for rs, (vc, vp) in LEAD_POS.items():
    hit = [r for r in c1 if r["vchr"] == vc and r["vpos"] == vp]
    nearest = min(c1, key=lambda r: abs(r["vpos"] - vp)) if c1 else None
    lead_query.append({
        "lead": rs, "pos": vp,
        "present_as_independent_signal": bool(hit),
        "subsets_with_hit": sorted({r["subset"] for r in hit}),
        "nearest_C1GALT1_signal": (
            {"subset": nearest["subset"], "variant": nearest["variant_id"],
             "pos": nearest["vpos"], "distance_bp": abs(nearest["vpos"] - vp)}
            if nearest else None),
    })

# ---- eGene 广度：每个轴基因在 28 个亚型中出现的次数 ----
breadth = {}
for gene in AXIS.values():
    subs = sorted({s for (s, g) in rows if g == gene})
    breadth[gene] = {"n_subsets_eGene": len(subs), "subsets": subs}

# ---- B 系亚型结论 ----
b_conc = {}
for gene in AXIS.values():
    b_conc[gene] = {
        s: next((o for o in out if o["cell_subset"] == s and o["gene"] == gene), None)
        for s in B_SUBSETS
    }

# ---- 写出 ----
tsv = os.path.join(RES, "TableS12_EAS_eQTL_crossref.tsv")
cols = ["cell_subset", "gene", "is_eGene_FDR0.05", "n_independent_signals",
        "top_variant_rank0", "top_variant_pos", "top_nominal_P", "top_slope",
        "min_nominal_P_any_signal", "B_lineage"]
with open(tsv, "w", encoding="utf-8") as f:
    f.write("\t".join(cols) + "\n")
    for o in out:
        f.write("\t".join(str(o[c]) for c in cols) + "\n")

summary = {
    "source": "ImmuNexUT E-GEAD-398 conditional_eQTL_FDR0.05 (Ota 2021 Cell; 416 Japanese; 28 subsets; hg38)",
    "n_subsets": 28,
    "n_axis_gene_rows": sum(subset_hits.values()),
    "axis_eGene_breadth": breadth,
    "B_lineage_conclusions": b_conc,
    "lead_variant_query": lead_query,
}
with open(os.path.join(RES, "step2_summary.json"), "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2, ensure_ascii=False)

print("=== 轴基因 eGene 广度（28 亚型中）===")
for g, d in breadth.items():
    print(f"  {g:12s} eGene in {d['n_subsets_eGene']:2d}/28 subsets")
print("\n=== B 系亚型（Naive_B/USM_B/SM_B/DN_B/Plasmablast）===")
for g in AXIS.values():
    for s in B_SUBSETS:
        o = next((x for x in out if x["cell_subset"] == s and x["gene"] == g), None)
        if o:
            print(f"  {g:12s} {s:12s} signals={o['n_independent_signals']:3d} "
                  f"top={o['top_variant_rank0']:<18s} pos={o['top_variant_pos']:<9d} "
                  f"P={o['top_nominal_P']:.3g} slope={o['top_slope']:+.3f}")
print("\n=== lead 位点查询（C1GALT1）===")
for q in lead_query:
    print(f"  {q['lead']} @{q['pos']}: present_as_independent_signal={q['present_as_independent_signal']}"
          f"  nearest={q['nearest_C1GALT1_signal']}")
print(f"\n[done] C1 step2 → {tsv}")
