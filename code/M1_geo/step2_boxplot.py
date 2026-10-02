# -*- coding: utf-8 -*-
"""
M1 · step2 箱线图 (5机制基因 × 3数据集, IgAN vs 对照)
数据来源: step1 输出的 {gse}_5genes_expr.csv (统一 log2/归一化尺度, 多探针均值)
显著性标注: 对应 *_best.csv 的 p 值 (Welch t-test, 多探针取最显著)
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
# 中文字体支持 (Windows: 微软雅黑/黑体)
_zh_ok = False
for _f in ("Microsoft YaHei", "SimHei", "SimSun"):
    try:
        font_manager.findfont(_f, fallback_to_default=False)
        plt.rcParams["font.sans-serif"] = [_f, "DejaVu Sans"]
        _zh_ok = True
        break
    except Exception:
        continue
plt.rcParams["axes.unicode_minus"] = False
if not _zh_ok:
    print("  [warn] 无中文字体, 图内中文将显示为方块")

OUT = _paths.at("阶段1/M1_GEO表达锚点")
FIG = OUT("figures")
os.makedirs(FIG, exist_ok=True)

GENES = ["C1GALT1", "C1GALT1C1", "ST6GALNAC2", "GALNT2", "GALNT12"]

# gse -> (标题, 分组tsv, 表达csv, best csv)
DATASETS = [
    ("GSE73953",  "PBMC · IgAN vs 健康对照(pooled)",      "GSE73953_groups.tsv",
     "GSE73953_5genes_expr.csv", "GSE73953_IgAN_vs_HC_best.csv"),
    ("GSE115857", "肾活检bulk · IgAN vs 活体供肾",        "GSE115857_groups.tsv",
     "GSE115857_5genes_expr.csv", "GSE115857_IgAN_vs_Ctrl_LD_best.csv"),
    ("GSE93798",  "肾小球 · IgAN vs 对照",                "GSE93798_groups.tsv",
     "GSE93798_5genes_expr.csv", "GSE93798_IgAN_vs_Ctrl_best.csv"),
]

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

def load_res(gse, grp_rel, expr_rel, best_rel):
    p_grp  = OUT(grp_rel)
    p_expr = OUT(expr_rel)
    p_best = OUT(best_rel)
    if not all(os.path.exists(x) for x in (p_grp, p_expr, p_best)):
        print(f"  [skip] 缺少 {gse} 输入文件"); return None
    grp  = pd.read_csv(p_grp, sep="\t", index_col=0, header=None, names=["grp"])
    expr = pd.read_csv(p_expr, index_col=0)
    best = pd.read_csv(p_best)
    return grp, expr, best

def plot_all(gse, title, case_name, ctrl_name, grp, expr, best):
    fig, axes = plt.subplots(1, 5, figsize=(21, 4.6))
    fig.suptitle(f"{gse} — {title}", fontsize=13, y=1.0)
    pval = best.set_index("gene")["pvalue"]
    logfc = best.set_index("gene")["logFC"]
    for ax, gene in zip(axes, GENES):
        if gene not in expr.index:
            ax.set_title(f"{gene}\n(no probe)", color="gray", fontsize=10)
            ax.axis("off"); continue
        s = expr.loc[gene]
        ca = np.asarray([s[c] for c in expr.columns if grp.loc[c, "grp"] == case_name], float)
        co = np.asarray([s[c] for c in expr.columns if grp.loc[c, "grp"] == ctrl_name], float)
        if len(ca) < 1 or len(co) < 1:
            ax.set_title(f"{gene}\n(缺组)", color="gray", fontsize=10); ax.axis("off"); continue
        bp = ax.boxplot([ca, co], tick_labels=[f"{case_name}\nn={len(ca)}", f"{ctrl_name}\nn={len(co)}"],
                        widths=0.5, patch_artist=True, showfliers=False,
                        medianprops=dict(color="black", lw=1.4),
                        boxprops=dict(lw=0.8), whiskerprops=dict(lw=0.8),
                        capprops=dict(lw=0.8))
        for patch, col in zip(bp["boxes"], ["#d62728", "#1f77b4"]):
            patch.set_facecolor(col); patch.set_alpha(0.72)
        for x, v in [(1, ca), (2, co)]:
            ax.scatter(np.random.default_rng(7).normal(x, 0.06, len(v)), v,
                       s=12, color="black", alpha=0.4, zorder=3, linewidths=0)
        # p 值标注
        if gene in pval.index:
            p = pval.loc[gene]
            ymax = max(np.nanmax(ca), np.nanmax(co)) if len(ca) and len(co) else np.nan
            ymin = min(np.nanmin(ca), np.nanmin(co))
            span = (ymax - ymin) if ymax != ymin else 1.0
            ax.text(0.5, 1.02, f"{fmt_p(p)} {star(p)}", transform=ax.transAxes,
                    ha="center", va="bottom", fontsize=9,
                    color="#c0392b" if p < 0.05 else "gray")
            fc = logfc.loc[gene]
            ax.set_title(f"{gene}\nlogFC={fc:+.2f}", fontsize=10)
        ax.set_ylabel("expression", fontsize=8)
        ax.tick_params(axis="x", labelsize=8)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, f"boxplot_{gse}.{ext}"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved figures/boxplot_{gse}.png/.pdf")

def main():
    for gse, title, grp_rel, expr_rel, best_rel in DATASETS:
        print("==", gse)
        if gse == "GSE73953":
            case_name, ctrl_name = "IgAN", "HC"
        elif gse == "GSE115857":
            case_name, ctrl_name = "IgAN", "Ctrl_LD"
        else:
            case_name, ctrl_name = "IgAN", "Ctrl"
        res = load_res(gse, grp_rel, expr_rel, best_rel)
        if res is None:
            continue
        grp, expr, best = res
        plot_all(gse, title, case_name, ctrl_name, grp, expr, best)

if __name__ == "__main__":
    main()
