# -*- coding: utf-8 -*-
"""
修复版 Figure_S1 (GSE73953/115857/93798 三数据集 × 5 基因 boxplot)
复用 step2_boxplot.py 主体，输出到复现包目录。
v3.2 修正：把 P 值与星号合并进 set_title（避免 ax.text 被 set_title 覆盖导致 P 值缺失）
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
from matplotlib import font_manager
for _f in ("Microsoft YaHei", "SimHei", "SimSun"):
    try:
        font_manager.findfont(_f, fallback_to_default=False)
        plt.rcParams["font.sans-serif"] = [_f, "DejaVu Sans"]
        break
    except Exception:
        continue
plt.rcParams["axes.unicode_minus"] = False

OUT = _paths.at("阶段1/M1_GEO表达锚点")
DST = _paths.at("阶段4.5_复现包/submission/additional_file_1")
GENES = ["C1GALT1", "C1GALT1C1", "ST6GALNAC2", "GALNT2", "GALNT12"]


def fmt_p(p):
    if pd.isna(p):
        return "NA"
    if p < 1e-4:
        return f"p={p:.1e}"
    if p < 0.05:
        return f"p={p:.3f}"
    return f"p={p:.2f}"


def star(p):
    if pd.isna(p):
        return ""
    if p < 1e-4:
        return "****"
    if p < 1e-3:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


# gse -> (标题, 分组tsv, 表达csv, best csv, case, ctrl)
DATASETS = [
    ("GSE73953",  "PBMC: IgAN vs pooled healthy controls",
     "GSE73953_groups.tsv", "GSE73953_5genes_expr.csv",
     "GSE73953_IgAN_vs_HC_best.csv", "IgAN", "HC"),
    ("GSE115857", "kidney biopsy bulk: IgAN vs living-donor controls",
     "GSE115857_groups.tsv", "GSE115857_5genes_expr.csv",
     "GSE115857_IgAN_vs_Ctrl_LD_best.csv", "IgAN", "Ctrl_LD"),
    ("GSE93798",  "glomeruli: IgAN vs controls",
     "GSE93798_groups.tsv", "GSE93798_5genes_expr.csv",
     "GSE93798_IgAN_vs_Ctrl_best.csv", "IgAN", "Ctrl"),
]


def plot_one(gse, title, case_name, ctrl_name, grp, expr, best):
    fig, axes = plt.subplots(1, 5, figsize=(21, 4.6))
    fig.suptitle(f"{gse} — {title}", fontsize=13, y=1.0)
    pval = best.set_index("gene")["pvalue"]
    logfc = best.set_index("gene")["logFC"]
    for ax, gene in zip(axes, GENES):
        if gene not in expr.index:
            ax.set_title(f"$\\it{{{gene.replace('1','1')}}}$\n(no probe)", color="gray", fontsize=10)
            ax.axis("off")
            continue
        s = expr.loc[gene]
        ca = np.asarray([s[c] for c in expr.columns if grp.loc[c, "grp"] == case_name], float)
        co = np.asarray([s[c] for c in expr.columns if grp.loc[c, "grp"] == ctrl_name], float)
        if len(ca) < 1 or len(co) < 1:
            ax.set_title(f"$\\it{{{gene}}}$\n(no group)", color="gray", fontsize=10)
            ax.axis("off")
            continue
        bp = ax.boxplot(
            [ca, co],
            tick_labels=[f"{case_name}\nn={len(ca)}", f"{ctrl_name}\nn={len(co)}"],
            widths=0.5, patch_artist=True, showfliers=False,
            medianprops=dict(color="black", lw=1.4),
            boxprops=dict(lw=0.8), whiskerprops=dict(lw=0.8),
            capprops=dict(lw=0.8),
        )
        for patch, col in zip(bp["boxes"], ["#d62728", "#1f77b4"]):
            patch.set_facecolor(col)
            patch.set_alpha(0.72)
        for x, v in [(1, ca), (2, co)]:
            ax.scatter(
                np.random.default_rng(7).normal(x, 0.06, len(v)), v,
                s=12, color="black", alpha=0.4, zorder=3, linewidths=0,
            )
        if gene in pval.index:
            p = pval.loc[gene]
            fc = logfc.loc[gene]
            # v3.2 修正：把 P 值与星号合并进 set_title，避免 ax.text 被覆盖导致 P 值缺失
            # 星号规则与 Fig1b 图注一致：仅当 |logFC|>1 且 P<0.05 才显示星号
            sig = abs(fc) > 1 and p < 0.05
            ax.set_title(
                f"$\\it{{{gene}}}$\nlogFC={fc:+.2f},  {fmt_p(p)} {star(p) if sig else 'ns'}",
                fontsize=9,
                color="black" if p >= 0.05 else "#c0392b",
            )
        ax.set_ylabel("expression", fontsize=8)
        ax.tick_params(axis="x", labelsize=8)
    fig.tight_layout()
    stem = DST(f"Figure_S1_{gse}_5genes_boxplot")
    fig.savefig(f"{stem}.png", dpi=200, bbox_inches="tight")
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {stem}.png/.pdf")


def main():
    for gse, title, grp_rel, expr_rel, best_rel, case, ctrl in DATASETS:
        print("==", gse)
        p_grp = OUT(grp_rel)
        p_expr = OUT(expr_rel)
        p_best = OUT(best_rel)
        if not all(os.path.exists(x) for x in (p_grp, p_expr, p_best)):
            print(f"  [skip] 缺少 {gse} 输入文件")
            continue
        grp = pd.read_csv(p_grp, sep="\t", index_col=0, header=None, names=["grp"])
        expr = pd.read_csv(p_expr, index_col=0)
        best = pd.read_csv(p_best)
        plot_one(gse, title, case, ctrl, grp, expr, best)


if __name__ == "__main__":
    main()