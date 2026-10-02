# -*- coding: utf-8 -*-
"""M2a 可视化：5基因 × 5语境 eQTL 效应与显著性
Fig A: 点图（每基因一横条语境，beta 定位，显著着色）
Fig B: -log10(p_beta) 热图（no-hit 置灰）
输入: <REPO>/data/M2a_eqtl/eQTL_5genes_x_context.tsv
输出: <REPO>/data/M2a_eqtl/figures/eQTL_5genes_forest.png/pdf
      <REPO>/data/M2a_eqtl/figures/eQTL_5genes_heatmap.png/pdf
"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import os, csv, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT_DIR = _paths.at("阶段1/M2a_eQTL语境")        # → <REPO>/data/M2a_eqtl
FIG_DIR = OUT_DIR("figures")
os.makedirs(FIG_DIR, exist_ok=True)

GENES = ["C1GALT1", "C1GALT1C1", "GALNT2", "GALNT12", "ST6GALNAC2"]
CTX = ["OneK1K_B_naive", "OneK1K_B_memory", "OneK1K_B_intermediate",
       "GTEx_v10_blood", "GTEx_v10_kidney_cortex"]
CTX_LABEL = {"OneK1K_B_naive": "OneK1K\nB_naive", "OneK1K_B_memory": "OneK1K\nB_memory",
             "OneK1K_B_intermediate": "OneK1K\nB_intermediate",
             "GTEx_v10_blood": "GTEx v10\nWhole Blood", "GTEx_v10_kidney_cortex": "GTEx v10\nKidney Cortex"}

def load():
    d = {}
    with open(OUT_DIR("eQTL_5genes_x_context.tsv"), encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            d[(r["gene"], r["dataset"])] = r
    return d

def pf(x):
    try:
        return float(x)
    except Exception:
        return None

data = load()
# ---------- Fig A: 点图 ----------
fig, axes = plt.subplots(len(GENES), 1, figsize=(9, 2.4 * len(GENES)), sharex=True)
for ax, gene in zip(axes, GENES):
    xs, ys, sigs, labels = [], [], [], []
    for i, ct in enumerate(CTX):
        r = data.get((gene, ct))
        b = pf(r["beta"]) if r and r["beta"] not in ("", "NA") else None
        p = pf(r["p_beta"]) if r and r["p_beta"] not in ("", "NA") else None
        ys.append(len(CTX) - i)  # 从上到下 B_naive..kidney
        if b is None:
            xs.append(0); sigs.append("none"); labels.append(CTX_LABEL[ct])
        else:
            xs.append(b)
            sigs.append("sig" if (p is not None and p < 0.05) else "ns")
            labels.append(CTX_LABEL[ct])
    for x, y, s, lab in zip(xs, ys, sigs, labels):
        if s == "none":
            ax.scatter(x, y, marker="x", s=60, color="gray", linewidths=1.2, label="no eQTL record")
        elif s == "sig":
            ax.scatter(x, y, marker="o", s=120, color="#C0392B" if x > 0 else "#1F618D", zorder=3)
            ax.text(x, y + 0.16, f"{x:+.2f}", ha="center", fontsize=8, color="#333")
        else:
            ax.scatter(x, y, marker="o", s=70, facecolors="none", edgecolors="#999", linewidths=1.0)
    ax.axvline(0, color="#bbb", lw=0.8, ls="--")
    ax.set_yticks(range(1, len(CTX) + 1))
    ax.set_yticklabels([CTX_LABEL[c] for c in CTX], fontsize=8)
    ax.set_ylim(0.3, len(CTX) + 0.8)
    ax.set_title(gene, loc="left", fontsize=11, fontweight="bold")
    ax.set_xlim(-1.2, 1.2)
axes[-1].set_xlabel("cis-eQTL effect (beta, normalized expression)", fontsize=9)
axes[0].legend(loc="upper right", fontsize=7, frameon=False)
fig.suptitle("M2a  eQTL effect of O-glycosylation axis genes across B-cell / tissue contexts\n(eQTL Catalogue r8; red/blue solid = gene-level p<0.05)", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(os.path.join(FIG_DIR, "eQTL_5genes_forest.png"), dpi=300)
fig.savefig(os.path.join(FIG_DIR, "eQTL_5genes_forest.pdf"))
plt.close(fig)

# ---------- Fig B: 热图 -log10(p_beta) ----------
mat = np.zeros((len(GENES), len(CTX))); annot = np.empty((len(GENES), len(CTX)), dtype=object)
for gi, gene in enumerate(GENES):
    for ci, ct in enumerate(CTX):
        r = data.get((gene, ct))
        p = pf(r["p_beta"]) if r and r["p_beta"] not in ("", "NA") else None
        b = pf(r["beta"]) if r and r["beta"] not in ("", "NA") else None
        if p is None:
            mat[gi, ci] = np.nan
            annot[gi, ci] = "no hit"
        else:
            mat[gi, ci] = min(-math.log10(max(p, 1e-300)), 12)
            annot[gi, ci] = f"{b:+.2f}\np={p:.1e}" if abs(p) < 0.05 else f"{b:+.2f}"
fig, ax = plt.subplots(figsize=(9, 4.2))
cmap = matplotlib.colormaps["Reds"].copy(); cmap.set_bad("#EFEFEF")
im = ax.imshow(mat, cmap=cmap, aspect="auto", vmin=0, vmax=8)
ax.set_xticks(range(len(CTX))); ax.set_xticklabels([CTX_LABEL[c] for c in CTX], fontsize=8)
ax.set_yticks(range(len(GENES))); ax.set_yticklabels(GENES, fontsize=10)
for gi in range(len(GENES)):
    for ci in range(len(CTX)):
        t = annot[gi, ci]
        ax.text(ci, gi, t, ha="center", va="center", fontsize=7.5,
                color="white" if (not np.isnan(mat[gi, ci]) and mat[gi, ci] > 5) else "#333")
ax.set_title("cis-eQTL significance [−log10(p_beta)] for O-glycosylation genes\n(red depth = stronger eQTL; grey = no record in dataset)", fontsize=10)
cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
cb.set_label("-log10 p_beta (cap 8)", fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "eQTL_5genes_heatmap.png"), dpi=300)
fig.savefig(os.path.join(FIG_DIR, "eQTL_5genes_heatmap.pdf"))
plt.close(fig)
print("figures done:", os.listdir(FIG_DIR))
