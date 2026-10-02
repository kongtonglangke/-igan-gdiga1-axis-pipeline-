# -*- coding: utf-8 -*-
"""阶段2 主证据1 —— 可视化：lead Gd-IgA1 位点 × 细胞型 eQTL 归属
图1 (eQTL_lead_matrix_heat.png/pdf)：4 lead × 5 语境，-log10(axis_gene_min_p) 热条
    底色规则：LEAD_ABSENT=深灰 / 未检测(gene not tested)=浅灰斜纹 / 其余按 -log10 p 渐变
图2 (C1GALT1_leads_beta.png/pdf)：C1GALT1 两个 lead 在各语境的 beta 点图（效应方向+CI 近似）
数据：lead_eQTL_matrix.tsv（步骤B 输出）
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
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = _paths.at("阶段2/主证据1_lead_eQTL")
FIG = HERE("figures")
os.makedirs(FIG, exist_ok=True)

LEADS = ["rs13226913", "rs10238682", "rs7856182", "rs5910940"]
LEAD_LABEL = {"rs13226913": "rs13226913\n(C1GALT1, Kiryluk'17)",
              "rs10238682": "rs10238682\n(C1GALT1, Wang'21)",
              "rs7856182": "rs7856182\n(GALNT12 下游19kb, Wang'21)",
              "rs5910940": "rs5910940\n(C1GALT1C1 基因内, Kiryluk'17)"}
DSET_ORDER = ["OneK1K_B_naive", "OneK1K_B_memory", "OneK1K_B_intermediate",
              "GTEx_v10_blood", "GTEx_v10_kidney_cortex"]
DSET_LABEL = {"OneK1K_B_naive": "OneK1K B naive\n(n=876)",
              "OneK1K_B_memory": "OneK1K B memory\n(n=726)",
              "OneK1K_B_intermediate": "OneK1K B inter.\n(n=618)",
              "GTEx_v10_blood": "GTEx v10 whole blood\n(n=853)",
              "GTEx_v10_kidney_cortex": "GTEx v10 kidney cortex\n(n=106)"}

rows = [l.rstrip("\n").split("\t") for l in open(HERE("lead_eQTL_matrix.tsv"), encoding="utf-8")]
rows = [r for r in rows if r and not r[0].startswith("lead_rsid")]
idx = {r[0]: i for i, r in enumerate(rows)}

def get(lead, ds):
    for r in rows:
        if r[0] == lead and r[5] == ds:
            return r
    return None

# ============ 图1：显著性热条 ============
fig, ax = plt.subplots(figsize=(9.5, 4.6))
import matplotlib.patches as mpatches
X = np.arange(len(DSET_ORDER))
width = 0.19
colors_ds = {"OneK1K_B_naive": "#2166ac", "OneK1K_B_memory": "#4393c3",
             "OneK1K_B_intermediate": "#92c5de", "GTEx_v10_blood": "#d6604d",
             "GTEx_v10_kidney_cortex": "#b2182b"}
handles = []
for li, lead in enumerate(LEADS):
    xpos = X + (li - 1.5) * width
    for xi, ds in enumerate(DSET_ORDER):
        r = get(lead, ds)
        axv = xpos[xi]
        if r is None:
            continue
        ap = r[9]
        if ap == "LEAD_ABSENT":
            ax.bar(axv, 1.0, width, color="#969696", edgecolor="k", lw=0.4)
            ax.text(axv, 0.5, "SNP\nabsent", ha="center", va="center", fontsize=6.5, color="white")
        elif ap == "":
            ax.bar(axv, 1.0, width, color="#d9d9d9", edgecolor="k", lw=0.4, hatch="//")
            ax.text(axv, 0.5, "gene\nnot tested", ha="center", va="center", fontsize=6.5, color="#666666")
        else:
            p = float(ap)
            lp = min(-np.log10(p), 12)
            ax.bar(axv, lp, width, color=colors_ds[ds], edgecolor="k", lw=0.4)
            sig = "***" if p < 5e-8 else ("**" if p < 1e-4 else ("*" if p < 0.05 else ""))
            if sig:
                ax.text(axv, lp + 0.25, sig, ha="center", va="bottom", fontsize=9, fontweight="bold")
    # legend 用 ds 颜色
for ds in DSET_ORDER:
    handles.append(mpatches.Patch(color=colors_ds[ds], label=DSET_LABEL[ds]))
handles += [mpatches.Patch(color="#969696", label="lead SNP 不在该数据集"),
            mpatches.Patch(color="#d9d9d9", hatch="//", label="注释基因未被检测 (B 细胞低表达)")]
ax.axhline(-np.log10(5e-8), color="red", ls="--", lw=0.9)
ax.text(len(DSET_ORDER) - 0.2, -np.log10(5e-8) + 0.25, "genome-wide sig (5e-8)", color="red", fontsize=7.5, ha="right")
ax.axhline(-np.log10(1e-4), color="gray", ls=":", lw=0.8)
ax.set_xticks(X)
ax.set_xticklabels([DSET_LABEL[d] for d in DSET_ORDER], fontsize=8.5)
ax.set_ylabel("-log10(p)  (lead SNP 对注释基因的最小名义 P)", fontsize=9.5)
ax.set_ylim(0, 13.5)
ax.set_title("Gd-IgA1 lead SNP 本尊 × 语境 cis-eQTL 归属\n(主证据1: B细胞/浆细胞机制归属)", fontsize=11.5)
ax.legend(handles=handles, loc="upper right", fontsize=7.5, ncol=2, frameon=True, framealpha=0.9)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "eQTL_lead_matrix_heat.png"), dpi=300)
fig.savefig(os.path.join(FIG, "eQTL_lead_matrix_heat.pdf"))
plt.close(fig)

# ============ 图2：C1GALT1 lead 的 beta 跨语境 ============
leads_c1 = ["rs13226913", "rs10238682"]
fig, ax = plt.subplots(figsize=(8.5, 3.8))
for li, lead in enumerate(leads_c1):
    off = -0.15 if li == 0 else 0.15
    ypos = []
    ylab = []
    betas, ses, ps = [], [], []
    for yi, ds in enumerate(DSET_ORDER):
        r = get(lead, ds)
        if r is None or r[9] == "" or r[9] == "LEAD_ABSENT":
            continue
        ap, ab, ase = r[9], r[10], r[11]
        if ap == "":
            continue
        p = float(ap); b = float(ab); s = float(ase) if ase else 0.1
        ypos.append(yi + off); betas.append(b); ses.append(max(s, 1e-6)); ps.append(p)
        ylab.append(DSET_ORDER[yi])
    ax.errorbar(betas, ypos, xerr=[1.96 * s for s in ses], fmt="o", ms=6,
                color="#2166ac" if li == 0 else "#b2182b", ecolor="0.3",
                elinewidth=1.2, capsize=3, label=LEAD_LABEL[lead].split("\n")[0] + " " + ("β: Kiryluk'17" if li==0 else "β: Wang'21"))
ax.axvline(0, color="k", lw=0.8)
ax.set_yticks(range(len(DSET_ORDER)))
ax.set_yticklabels([DSET_LABEL[d] for d in DSET_ORDER], fontsize=8.5)
ax.set_xlabel("beta (eQTL effect, per-allele)", fontsize=10)
ax.set_title("C1GALT1 lead SNP → C1GALT1 表达效应量跨语境对比\n(rs13226913 与 rs10238682 同基因体, 效应方向一致为 +)", fontsize=11)
ax.legend(fontsize=8, loc="lower right")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "C1GALT1_leads_beta.png"), dpi=300)
fig.savefig(os.path.join(FIG, "C1GALT1_leads_beta.pdf"))
plt.close(fig)

print("[OK] figures saved to", FIG)
