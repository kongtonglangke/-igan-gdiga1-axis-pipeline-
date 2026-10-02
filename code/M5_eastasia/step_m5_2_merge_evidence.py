#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""M5-2: 合并多源 → M5 东亚汇总中间表
源1: GCST90018866 EAF/beta/P (step_m5_1 输出)
源2: Wang 2021 汉族 Gd-IgA1 P (lead_x_Wang2021_GdIgA1_P.tsv)
源3: OneK1K B_naive MAF (lead_hits_all.tsv, 新加坡多族裔)
源4: 阶段2 主要 eQTL 语境轴_minP (lead_eQTL_matrix.tsv)
输出: M5_东亚/lead_x_东亚证据汇总.tsv
"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import csv, os

D = _paths.LEGACY
def rd(p, **kw):
    return list(csv.DictReader(open(p, encoding="utf-8"), delimiter="\t", **kw))

# 源1
ig = {r["lead_rsid"]: r for r in rd(D("阶段3/M5_东亚/lead_x_IgAN_eaf.tsv"))}
# 源2
wa = {r["lead_rsid"]: r["gda_P"] for r in rd(D("阶段2/主证据2_共享结构/lead_x_Wang2021_GdIgA1_P.tsv"))}
# 源3: OneK1K B_naive MAF（MAF 属 lead SNP 本身，取该 lead 在该语境的任意行）
hits = rd(D("阶段2/主证据1_lead_eQTL/lead_hits_all.tsv"))
ok = {}
for r in hits:
    if r["dataset"] == "OneK1K_B_naive":
        ok.setdefault(r["lead_rsid"], r["maf"])

rows = []
order = ["rs13226913", "rs10238682", "rs7856182", "rs5910940"]
for rs in order:
    g = ig.get(rs, {})
    rows.append({
        "lead_rsid": rs,
        "lead_gene": {"rs13226913": "C1GALT1", "rs10238682": "C1GALT1",
                      "rs7856182": "GALNT12", "rs5910940": "C1GALT1C1"}[rs],
        "igan_meta_EAF": g.get("EAF", "NA"),
        "igan_meta_beta": g.get("beta", "NA"),
        "igan_meta_P": g.get("p", "NA"),
        "Wang2021_Han_GdIgA1_P": wa.get(rs, "NA"),
        "OneK1K_Bnaive_MAF": ok.get(rs, "NA"),
    })

out = D("阶段3/M5_东亚/lead_x_东亚证据汇总.tsv")
with open(out, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows)
print("saved:", out)
for r in rows: print(r)
