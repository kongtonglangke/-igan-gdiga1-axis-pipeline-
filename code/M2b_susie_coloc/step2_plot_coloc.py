# -*- coding: utf-8 -*-
"""阶段2 主证据2 —— 可视化（主图 2-3 草稿）
图1 coloc_PPH_stacked.png/pdf：5 基因 × 5 语境 coloc.abf 五假设后验概率堆叠条
    直观展示: 可算组全部 PPH4≈0; PPH0 或 PPH1 占优(无 IgAN 侧共享信号)
图2 lead_phenotype_ladder.png/pdf：4 lead 在 Gd-IgA1(Wang'21, P-only) vs IgAN(Sakaue'21)
    的分层证据: 上半段 -log10(Gd-IgA1 P) 极高, 下半段 IgAN P>0.1 → 机制轴效应
    特异于 Gd-IgA1 产生、不传导至 IgAN 疾病易感。
数据: coloc_5genes_x_contexts.tsv, lead_x_Wang2021_GdIgA1_P.tsv, lead_x_IgAN_Pvalue.tsv
"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = _paths.at("阶段2/主证据2_共享结构")
FIG = HERE("figures")
os.makedirs(FIG, exist_ok=True)

GENE_ORDER = ["C1GALT1", "GALNT2", "GALNT12", "ST6GALNAC2", "C1GALT1C1"]
CTX_ORDER = ["GTEx_v10_blood", "GTEx_v10_kidney_cortex",
             "OneK1K_B_naive", "OneK1K_B_memory", "OneK1K_B_intermediate"]
CTX_LABEL = {"GTEx_v10_blood": "GTEx blood", "GTEx_v10_kidney_cortex": "GTEx kidney",
             "OneK1K_B_naive": "OneK1K B naive", "OneK1K_B_memory": "OneK1K B memory",
             "OneK1K_B_intermediate": "OneK1K B inter."}

# ============ 图1: coloc PPH 堆叠 ============
co = pd.read_csv(HERE("coloc_5genes_x_contexts.tsv"), sep="\t")
# 只取有 PPH 的组（n_SNP>=5 且 PPH0 非 NaN）
co = co[co["PPH0"].notna()].copy()
co["label"] = co["gene"] + " / " + co["context"].map(CTX_LABEL)

H_colors = {"PPH0": "#d9d9d9", "PPH1": "#4292c6", "PPH2": "#9ecae1",
            "PPH3": "#fdae6b", "PPH4": "#e31a1c"}
labels_all = [f"{g} / {CTX_LABEL[c]}" for g in GENE_ORDER for c in CTX_ORDER]
labels_ok = [l for l in labels_all if l in set(co["label"])]

fig, ax = plt.subplots(figsize=(10.5, 6.4))
left = np.zeros(len(labels_ok))
for h in ["PPH0", "PPH1", "PPH2", "PPH3", "PPH4"]:
    vals = co.set_index("label").loc[labels_ok, h].values
    ax.barh(range(len(labels_ok)), vals, left=left, color=H_colors[h],
            edgecolor="white", lw=0.3, label={"PPH0":"H0 no assoc","PPH1":"H1 eQTL only",
            "PPH2":"H2 IgAN only","PPH3":"H3 distinct","PPH4":"H4 shared (coloc)"}[h])
    left = left + vals
for i, l in enumerate(labels_ok):
    row = co.set_index("label").loc[l]
    n = int(row["n_SNP"])
    ax.text(1.012, i, f"n={n}", va="center", fontsize=7, color="#333333")
ax.set_yticks(range(len(labels_ok)))
ax.set_yticklabels(labels_ok, fontsize=8.5)
ax.set_xlim(0, 1.10)
ax.set_xlabel("Posterior probability of H0-H4 (coloc.abf, p1=p2=1e-4, p12=1e-5)", fontsize=9.5)
ax.set_title("eQTL (mechanism axis) vs IgAN GWAS (GCST90018866): colocalization PPH\n"
             "all computable loci: PPH4<0.031  →  no shared causal variant with IgAN susceptibility",
             fontsize=10.5)
ax.invert_yaxis()
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.09), ncol=5, fontsize=8.2,
          frameon=False)
fig.subplots_adjust(left=0.22, bottom=0.18, top=0.90)
fig.savefig(os.path.join(FIG, "coloc_PPH_stacked.png"), dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(FIG, "coloc_PPH_stacked.pdf"), bbox_inches="tight")
plt.close(fig)
print("saved coloc_PPH_stacked")

# ============ 图2: lead 分层证据 (Gd-IgA1 vs IgAN) ============
gda = pd.read_csv(HERE("lead_x_Wang2021_GdIgA1_P.tsv"), sep="\t")
igan = pd.read_csv(HERE("lead_x_IgAN_Pvalue.tsv"), sep="\t")
m = gda.merge(igan, on="lead_rsid")
m["gene"] = m["lead_gene_x"]
m["rs_gene"] = m["lead_rsid"] + "  (" + m["gene"] + ")"
order = ["rs13226913  (C1GALT1)", "rs10238682  (C1GALT1)",
         "rs7856182  (GALNT12)", "rs5910940  (C1GALT1C1)"]
m["rs_gene"] = pd.Categorical(m["rs_gene"], categories=order, ordered=True)
m = m.sort_values("rs_gene")
y = np.arange(len(m))

fig, ax = plt.subplots(figsize=(8.2, 4.6))
# Gd-IgA1: -log10 P (P-only)
lp_gda = -np.log10(m["gda_P"].clip(lower=1e-300))
# IgAN: P -> -log10; 未命中(X) 特殊处理
lp_igan = np.where(m["igan_p"].notna(), -np.log10(m["igan_p"].clip(lower=1e-300)), np.nan)
w = 0.36
b1 = ax.barh(y + w/2, lp_gda, height=w, color="#b2182b", alpha=0.9, label="Gd-IgA1  (Wang'21, P-only)")
b2 = ax.barh(y - w/2, lp_igan, height=w, color="#4393c3", alpha=0.9, label="IgAN disease  (Sakaue'21)")
# 未命中标注
for i in range(len(m)):
    if np.isnan(lp_igan[i]):
        ax.text(0.15, i - w/2, "chrX not in GWAS", va="center", fontsize=7, color="#666666")
ax.axvline(-np.log10(0.05), color="#666666", ls=":", lw=0.8)
ax.axvline(-np.log10(5e-8), color="#e31a1c", ls="--", lw=0.9)
ax.text(-np.log10(0.05)+0.05, len(m)-0.35, "P=0.05", fontsize=7, color="#666666")
ax.text(-np.log10(5e-8)+0.05, len(m)-0.35, "P=5e-8", fontsize=7, color="#e31a1c")
# P 值文本
for i in range(len(m)):
    if pd.notna(m["gda_P"].iloc[i]):
        ax.text(lp_gda[i]+0.12, i + w/2, f"{m['gda_P'].iloc[i]:.1e}", va="center", fontsize=6.8, color="#b2182b")
    if pd.notna(m["igan_p"].iloc[i]):
        ax.text(lp_igan[i]+0.12, i - w/2, f"{m['igan_p'].iloc[i]:.3f}", va="center", fontsize=6.8, color="#4393c3")
ax.set_yticks(y)
ax.set_yticklabels(m["rs_gene"], fontsize=8.8)
ax.set_xlabel("-log10(P)", fontsize=10)
ax.set_title("Published Gd-IgA1 lead SNPs: effect restricted to Gd-IgA1 level,\n"
             "NOT to IgAN disease susceptibility (mechanism axis = modifier, not driver)",
             fontsize=9.5, pad=14)
ax.legend(loc="lower right", fontsize=8)
ax.set_xlim(0, max(lp_gda.max()*1.25, 12))
fig.subplots_adjust(left=0.21, right=0.97, top=0.82, bottom=0.14)
fig.savefig(os.path.join(FIG, "lead_phenotype_ladder.png"), dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(FIG, "lead_phenotype_ladder.pdf"), bbox_inches="tight")
plt.close(fig)
print("saved lead_phenotype_ladder")
