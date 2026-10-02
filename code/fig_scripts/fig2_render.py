"""fig2_render.py — Fig.2 渲染脚本（2026-09-10 v3.4 成熟序重排版）
作者：AI 王严课题组代办｜依赖：matplotlib/pandas/numpy（fig43 venv）

读：
  阶段1/M2a_eQTL语境/heatmap/eQTL_5genes_summary_heat.tsv  (panel a 矩阵)
  阶段2/主证据1_lead_eQTL/lead_eQTL_matrix.tsv              (panel b 矩阵)
  阶段2/主证据1_lead_eQTL/lead_hits_all.tsv                  (panel c β±SE 备查)

写：
  阶段4/主图/out/Fig2_v3.4.{png,pdf,tiff}（投稿版；3 面板垂直堆叠）

版本链：
  v01/v02  2026-09-06 常规样式微调（图注对齐 R2 标题、a 灰格标签单行化、
        a/b x 轴标签去前导、b P 值文本替代 **、图例外置、c legend lower left、
        ×/β/− unicode 恢复）
  v1.0  2026-09-06 步骤 4.3 冻结
  v2.0  2026-09-06 重叠修复（图例外置右侧、灰格 hatch ||| + 白底文字、
        P=9.4e-6 紧贴框顶）
  v3.0  2026-09-07 用户反馈"Fig2 b 图第一个柱子颜色混了"修复：
        (1) 组内柱步长 = 柱宽(0.18) → 柱与柱零间隙紧贴 → 步长 0.24、
            柱宽 0.15，柱间隙 0.09 unit 颜色边界清晰；
        (2) rs13226913 B_naive 高亮框默认 vermillion = rs10238682 柱色 =
            GW 虚线色 → 蓝柱被朱红框包着与朱红柱紧贴"混色" →
            高亮框改黑色 #000000，GW 虚线独占朱红语义；
        (3) panel a 灰格文字 "B low expr." → "low expr."（10 字符
            7pt ≈ 0.87 unit 超灰格 0.8 宽左右各 3px → 8 字符 0.70 unit
            containment 归零）。
  v3.1  2026-09-07 用户二轮反馈：
        (a) panel a colorbar 太靠右（figure.colorbar + constrained_layout
            留 ~20% 空白）→ make_axes_locatable.append_axes("right",
            size="3%", pad="2%") 紧贴热图；
        (b) panel b P=9.4e-6 文本溢出 xlim=-0.5 左界 → 左扩 xlim -0.85；
        (c) panel b 删除首柱 highlight_box（用户要求与其他柱统一无框）；
        (d) panel c 0.0 线被底部 spine 压住 → ylim 下界 0→-0.04、
            线 lw 0.7→0.9 独立可见。
  v3.2  2026-09-07 用户三轮反馈（panel a 最底行刻度标签重叠）：
        panel a x 轴两行式标签 "(n=876, OneK1K)" 第二行 15 字符
        ≈ 156px > cell 中心距 122px → B_naive/B_memory/B_inter./
        GTEx blood 标签两两重叠 6~34px（tick labels 不属 ax.texts，
        check_text 曾漏检）→ 标签去掉 ", OneK1K" 后缀与 panel b
        统一格式 "B_naive\n(n=876)"（7 字符 ≈ 55px 留足余量）；
        _fig2_check.py 新增 check_ticklabels（tick-tick bbox 两两
        求交）防回归。
        自检 _fig2_check.py：text/patch/containment/ticklabels
        全 0 overlap + tick-vs-text 交叉 0 + panel b 柱色像素采样
        全 OK（0 mismatch）。
  v3.3  2026-09-09 投稿版（L7 整合轮）：panel a 标题补第二行副题
        "(red depth = stronger eQTL; grey = no cis-eQTL record)"，
        与 ch8 图注口径一致。
  v3.4  2026-09-10 用户反馈（Fig 2c x 轴序）——**panel c 改为成熟序**：
        b_states 由 naive→memory→intermediate 改为
        **naive→intermediate→memory**（= B2/L1b 与 Fig. S7 使用的
        STAGE_ORDER 口径）。数据本身不变、仅列序改变；改序后折线呈
        "高位——骤降——部分回升"（rs13226913: 0.174→0.035→0.081），
        与 ch8 Fig 2c 图例"非单调、呈降—回升"的诚实表述严格一致，
        亦与 B2 异质性检验（Cochran's Q P≈0.046; meta-regression
        P≈0.018）的结论同向。此前 naive→memory→intermediate 的列序会
        在视觉上制造一条并不存在的单调下降趋势；panel c x 序 = 成熟序，
        与 ch8 图例及 B2 异质性检验同向。
  v3.6  2026-09-10 图文一致性修复（阶段4 ↔ 阶段4.5 全量对齐）：panel a
        标题去掉 "Fig. 2 — " 图号前缀——全套 7 张主图的 panel 标题均不带
        图号（Fig.1/3/4/5/6/GA 一律 "a/b/c/d 主题"），此处为排版遗留。
        未改动任何数据、版式几何或结论。
"""
from __future__ import annotations

# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from style_common import (  # noqa: E402
    OKABE_ITO, C_RS13226913, C_RS10238682, C_RS7856182, C_RS5910940,
    C_NOTTEST, C_GW_LINE, HEAT_CMAP,
    FS_PANEL_LABEL, FS_AXIS_LABEL, FS_TICK, FS_LEGEND, FS_INCELL,
    new_fig_three_panels, save_v01, save_final, hatch_not_tested, highlight_box, mm_to_in,
)
from mpl_toolkits.axes_grid1 import make_axes_locatable  # noqa: E402

PROJECT_ROOT = _paths.LEGACY
F2A_SUMMARY = PROJECT_ROOT("阶段1/M2a_eQTL语境/heatmap/eQTL_5genes_summary_heat.tsv")
F2B_MATRIX  = PROJECT_ROOT("阶段2/主证据1_lead_eQTL/lead_eQTL_matrix.tsv")
F2C_HITS    = PROJECT_ROOT("阶段2/主证据1_lead_eQTL/lead_hits_all.tsv")
OUT_DIR     = PROJECT_ROOT("阶段4/主图/out")
os.makedirs(OUT_DIR, exist_ok=True)


def _fmt_p(p):
    try:
        p = float(p)
    except Exception:  # noqa: BLE001
        return "NA"
    if p < 1e-4:
        return f"{p:.1e}".replace("e-0", "e-")
    if p < 0.001:
        return f"{p:.1e}"
    if p < 0.01:
        return f"{p:.4f}"
    return f"{p:.3f}"


def panel_a(ax):
    """5 基因 × 5 语境 cis-eQTL 显著性热图。"""
    df = pd.read_csv(F2A_SUMMARY, sep="\t")
    genes = ["C1GALT1", "C1GALT1C1", "GALNT2", "GALNT12", "ST6GALNAC2"]
    # v3.2：x 轴标签去掉 ", OneK1K" 后缀（15 字符 ≈ 156px 超 cell 中心距
    # 122px → 相邻刻度标签两两重叠 34px）→ 与 panel b 标签格式统一
    # "B_naive\n(n=876)"，单行最短 "(n=876)" 7 字符 ≈ 55px 留足余量
    contexts = [
        ("OneK1K_B_naive",         "B_naive\n(n=876)"),
        ("OneK1K_B_memory",        "B_memory\n(n=726)"),
        ("OneK1K_B_intermediate",  "B_inter.\n(n=618)"),
        ("GTEx_v10_blood",         "GTEx blood\n(n=853)"),
        ("GTEx_v10_kidney_cortex", "GTEx kidney\n(n=106)"),
    ]
    cap = 8.0
    mat = np.full((len(genes), len(contexts)), np.nan)
    annot = [["" for _ in contexts] for _ in genes]
    for i, g in enumerate(genes):
        for j, (ctx, _) in enumerate(contexts):
            row = df[(df["gene"] == g) & (df["dataset"] == ctx)]
            if row.empty:
                continue
            r = row.iloc[0]
            if r["call"] == "no-hit":
                annot[i][j] = "no-hit"
                continue
            p_beta = float(r["p_beta"])
            if not np.isfinite(p_beta) or p_beta <= 0:
                v = 0.0
            else:
                v = min(-np.log10(p_beta), cap)
            mat[i, j] = v
            beta = float(r["beta"])
            annot[i][j] = f"{beta:+.2f}\np={_fmt_p(p_beta)}"

    im = ax.imshow(np.where(np.isnan(mat), 0, mat), cmap=HEAT_CMAP, vmin=0, vmax=cap,
                   aspect="auto", origin="upper")

    # v02: 灰格 patch 缩小到不遮蔽相邻 cell 数值（hatch 留 0.85 宽）
    for i, g in enumerate(genes):
        for j, (ctx, _) in enumerate(contexts):
            row = df[(df["gene"] == g) & (df["dataset"] == ctx)]
            if row.empty:
                continue
            r = row.iloc[0]
            if r["call"] == "no-hit":
                reason = "X-linked" if g == "C1GALT1C1" and ctx.startswith("OneK1K") else "low expr."
                # 灰底 + 稀疏竖纹（"|" 不与水平文字笔画冲突）
                ax.add_patch(plt.Rectangle((j - 0.40, i - 0.40), 0.80, 0.80,
                                           facecolor="#F0F0F0",
                                           edgecolor=OKABE_ITO["gray"],
                                           hatch="|||", lw=0.4, zorder=1))
                # 文本放 cell 中央 + 半透明白底（彻底脱离斜纹干扰）
                ax.text(j, i, reason, ha="center", va="center",
                        fontsize=FS_INCELL - 0.5, color="black",
                        weight="bold", zorder=4,
                        bbox=dict(boxstyle="round,pad=0.2",
                                  facecolor="white", edgecolor=OKABE_ITO["gray"],
                                  lw=0.4, alpha=0.92))

    for i in range(len(genes)):
        for j in range(len(contexts)):
            if annot[i][j] == "no-hit":
                continue
            t = annot[i][j]
            color = "white" if (not np.isnan(mat[i, j]) and mat[i, j] >= 5) else "black"
            ax.text(j, i, t, ha="center", va="center",
                    fontsize=FS_INCELL, color=color)

    ax.set_xticks(range(len(contexts)))
    ax.set_xticklabels([c[1] for c in contexts], fontsize=FS_TICK)
    ax.set_yticks(range(len(genes)))
    genes_italic = [r"$\it{" + g + "}$" for g in genes]
    ax.set_yticklabels(genes_italic, fontsize=FS_TICK)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    ax.set_title("a   5 axis genes × 5 contexts (cis-eQTL)\n"
                 "(red depth = stronger eQTL; grey = no cis-eQTL record)",
                 fontsize=FS_PANEL_LABEL, loc="left", pad=4, weight="bold")
    # v3.1：colorbar 改 make_axes_locatable 紧贴热图右侧（之前 figure.colorbar
    # 配 constrained_layout 让色条与 axes 间留出 ~20% 宽空白，热图被推左
    # → "刻度表太靠右"），3% 宽 + 2% 间距，紧凑且与图例一致
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad="2%")
    cbar = ax.figure.colorbar(im, cax=cax)
    cbar.ax.tick_params(labelsize=FS_TICK - 0.5)
    cbar.set_label("-log10(p_beta)", fontsize=FS_TICK - 0.5)


def panel_b(ax):
    """4 lead × 5 语境 cis-eQTL 复现矩阵（−log10P 条形 + 不可测灰格）。

    v3.0（2026-09-07 用户反馈"第一个柱子颜色混了"）：
      (1) 组内柱中心步长 0.18(=柱宽)→0.24、柱宽 0.18→0.15：
          相邻柱零间隙紧贴 → 间隙 0.09 unit，颜色边界清晰分离；
      (2) highlight_box 默认色 C_GW_LINE(vermillion)=rs10238682 柱色 =
          GW 虚线色，蓝柱(rs13226913)被朱红框包着又与右侧朱红柱零距
          紧贴 → 视觉"第一个柱子颜色混了" → 高亮框改黑色 (#000000)，
          GW 虚线独占朱红语义，不再与任何柱撞色；
      (3) NA 灰格宽 0.18→0.15 与新柱宽一致。

    v3.1（2026-09-07 用户反馈）：
      (a) "P=9.4e-6" 文本在第一柱上方，文本宽 ~0.55 unit 但 xlim 左界
          -0.5 容不下 → 左扩 xlim 至 -0.85，确保 P 文本完全在 axes 内；
      (b) 首柱(蓝柱)外加黑色 highlight_box 让其与红框语义断裂，
          用户要求"统一没有" → 删除 highlight_box 调用，靠 P 文本
          + 黑色描边（所有柱都有 edgecolor="black" lw=0.4）自然标注。
    """
    df = pd.read_csv(F2B_MATRIX, sep="\t")
    leads = [
        ("rs13226913", "C1GALT1",   C_RS13226913),
        ("rs10238682", "C1GALT1",   C_RS10238682),
        ("rs7856182",  "GALNT12",   C_RS7856182),
        ("rs5910940",  "C1GALT1C1", C_RS5910940),
    ]
    contexts = [
        ("OneK1K_B_naive",         "B_naive\n(n=876)"),
        ("OneK1K_B_memory",        "B_memory\n(n=726)"),
        ("OneK1K_B_intermediate",  "B_inter.\n(n=618)"),
        ("GTEx_v10_blood",         "GTEx blood\n(n=853)"),
        ("GTEx_v10_kidney_cortex", "GTEx kidney\n(n=106)"),
    ]
    n_lead = len(leads); n_ctx = len(contexts)
    bar_w = 0.15   # 柱宽（v3.0: 0.18→0.15）
    step  = 0.24   # 组内步长（v3.0: 0.18→0.24，柱间隙 0.09）
    ax.axhline(-np.log10(5e-8), color=C_GW_LINE, lw=1.2, ls="--", zorder=2)
    # GW 线标签：右下角（避开所有柱；图例在右上外）
    ax.text(0.98, 0.96,
            "genome-wide sig  (5×10⁻⁸)", color=C_GW_LINE,
            fontsize=FS_TICK - 0.5, ha="right", va="top",
            transform=ax.transAxes,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                      edgecolor="none", alpha=0.92))

    for i, (rsid, gene, color) in enumerate(leads):
        for j, (ctx, _) in enumerate(contexts):
            row = df[(df["lead_rsid"] == rsid) & (df["dataset"] == ctx)]
            if row.empty:
                continue
            r = row.iloc[0]
            x = j + (i - 1.5) * step
            if str(r["axis_gene_tested"]).strip().lower() != "yes":
                ax.add_patch(plt.Rectangle((x - bar_w/2, 0), bar_w, 0.3,
                                           facecolor="#EEEEEE",
                                           edgecolor=C_NOTTEST,
                                           hatch="///", lw=0.6, zorder=1))
                ax.text(x, 0.18, "NA", ha="center", va="bottom",
                        fontsize=FS_INCELL - 1, color=C_NOTTEST)
                continue
            p = float(r["axis_min_p"])
            v = -np.log10(p) if p > 0 else 0.0
            ax.bar(x, v, width=bar_w, color=color, edgecolor="black", lw=0.4, zorder=3)

        # rs13226913 B_naive P 注释（v3.1 删除 highlight_box，仅靠 P 文本+柱自身
        # 黑色描边统一标注，与其他柱视觉一致）
        if rsid == "rs13226913":
            row = df[(df["lead_rsid"] == rsid) & (df["dataset"] == "OneK1K_B_naive")]
            if not row.empty:
                p = float(row.iloc[0]["axis_min_p"])
                x_star = 0 + (i - 1.5) * step
                v_star = -np.log10(p)
                # P 文本紧贴框上方 + 留白（远离 GW 线，柱顶 v=5.22 < GW 7.3）
                ax.text(x_star, v_star + 0.30, f"P={_fmt_p(p)}",
                        ha="center", va="bottom",
                        fontsize=FS_TICK, color="black", weight="bold")

    # 图例外置：右上 axes 外（bbox_inches='tight' 自动扩展画布）
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=c, edgecolor="black", lw=0.4)
               for _, _, c in leads]
    labels = [f"{rs} ({g})" for rs, g, _ in leads]
    ax.legend(handles, labels,
              loc="upper left", bbox_to_anchor=(1.02, 1.00),
              fontsize=FS_LEGEND,
              frameon=True, framealpha=0.95, ncol=1,
              title="lead × gene", title_fontsize=FS_LEGEND,
              borderaxespad=0.4)

    ax.set_xticks(range(n_ctx))
    ax.set_xticklabels([c[1] for c in contexts], fontsize=FS_TICK)
    # v3.1：xlim 左扩 -0.5→-0.85 给首柱上方 "P=9.4e-6" 文本（宽 ~0.55 unit，
    # 中心 x=-0.36）完整放入 axes 内
    ax.set_xlim(-0.85, n_ctx - 0.4)
    ax.set_ylim(0, 13.5)
    ax.set_ylabel(r"$-\log_{10}(P)$", fontsize=FS_AXIS_LABEL)
    ax.tick_params(axis="y", labelsize=FS_TICK)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    ax.set_title("b   4 Gd-IgA1 lead SNPs × 5 contexts — cis-eQTL replication",
                 fontsize=FS_PANEL_LABEL, loc="left", pad=4, weight="bold")


def panel_c(ax):
    """C1GALT1 leads β 沿 B 分化轨迹（rs13226913 & rs10238682; n=876/618/726）。

    v3.4：x 轴改为**成熟序** naive → intermediate → memory（此前为
    naive → memory → intermediate，会视觉制造单调下降的假趋势）。
    """
    df = pd.read_csv(F2C_HITS, sep="\t")
    b_states = [
        ("OneK1K_B_naive",         "B_naive (n=876)"),
        ("OneK1K_B_intermediate",  "B_intermediate (n=618)"),
        ("OneK1K_B_memory",        "B_memory (n=726)"),
    ]
    xs = list(range(len(b_states)))

    for rsid, color, label in [
        ("rs13226913", C_RS13226913, r"$\it{rs13226913}$ (Kiryluk'17)"),
        ("rs10238682", C_RS10238682, r"$\it{rs10238682}$ (Wang'21)"),
    ]:
        ys, ses = [], []
        for ctx, _ in b_states:
            sub = df[(df["lead_rsid"] == rsid) &
                     (df["dataset"] == ctx) &
                     (df["eQTL_gene_symbol"] == "C1GALT1")]
            if sub.empty:
                ys.append(np.nan); ses.append(np.nan); continue
            r = sub.iloc[0]
            ys.append(float(r["beta"]))
            ses.append(float(r["se"]))
        ys = np.array(ys); ses = np.array(ses)
        ax.errorbar(xs, ys, yerr=ses, fmt="o-", color=color, lw=1.5,
                    capsize=3, markersize=6, label=label, zorder=4)

    # v02: P 标签右上
    p_naive_rs13226913 = df[(df["lead_rsid"] == "rs13226913") &
                            (df["dataset"] == "OneK1K_B_naive") &
                            (df["eQTL_gene_symbol"] == "C1GALT1")].iloc[0]["pvalue"]
    p_naive_rs10238682 = df[(df["lead_rsid"] == "rs10238682") &
                            (df["dataset"] == "OneK1K_B_naive") &
                            (df["eQTL_gene_symbol"] == "C1GALT1")].iloc[0]["pvalue"]
    ax.text(0.96, 0.96,
            f"P={_fmt_p(p_naive_rs13226913)}  (Kiryluk'17, B_naive)", color=C_RS13226913,
            fontsize=FS_TICK, ha="right", va="top", weight="bold",
            transform=ax.transAxes)
    ax.text(0.96, 0.85,
            f"P={_fmt_p(p_naive_rs10238682)}  (Wang'21, B_naive)", color=C_RS10238682,
            fontsize=FS_TICK, ha="right", va="top",
            transform=ax.transAxes)

    # v3.1（2026-09-07 用户反馈）：0.0 线被底部 spine 压住 → ylim 下界
    # 由 (0,...) 改为 (-0.04,...)，让底部 spine 落到 y=-0.04、0.0 线在
    # y=0 处独立可见；同时 y-tick "0.0" 自然落到 0.0 线位置（图文一致）
    ax.axhline(0, color="black", lw=0.9, zorder=2)
    ax.set_xticks(xs)
    ax.set_xticklabels([b[1].split(" ")[0] + "\n" + b[1].split(" ", 1)[1] for b in b_states],
                       fontsize=FS_TICK)
    ax.set_xlim(-0.3, len(b_states) - 0.7)
    ax.set_ylim(-0.04, 0.25)
    ax.set_ylabel(r"$\beta$ (per-allele eQTL effect on $\it{C1GALT1}$)",
                 fontsize=FS_AXIS_LABEL)
    ax.tick_params(axis="y", labelsize=FS_TICK)
    # v02: legend 改 lower left
    ax.legend(loc="lower left", fontsize=FS_LEGEND, frameon=True, framealpha=0.95)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    ax.set_title("c   C1GALT1 lead eQTL $\\beta$ along B differentiation",
                 fontsize=FS_PANEL_LABEL, loc="left", pad=4, weight="bold")


def main():
    fig, axes = new_fig_three_panels(panel_h_mm=(60, 70, 60))
    panel_a(axes[0])
    panel_b(axes[1])
    panel_c(axes[2])
    # v2.0：图-文整合后第二轮重排——图例外置、灰格文本可读性、P 注释避 GW
    out_stem = os.path.join(OUT_DIR, os.environ.get("OUT_STEM", "Fig2_v3.4"))
    save_final(fig, out_stem)
    plt.close(fig)
    print("WROTE", out_stem + ".png/.pdf/.tiff")


if __name__ == "__main__":
    main()