"""fig5_render.py — Fig.5 东亚异质性（2026-09-06 步骤 4.3 制作第五图）
作者：AI 王严课题组代办｜依赖：matplotlib/pandas/numpy（fig43 venv）

读：
  阶段3/M5_东亚/M5_gnomAD_EASEUR_频率_20260903.tsv  (panel a 4 lead × 3 人群 AF)
  阶段3/M5_东亚/lead_x_东亚证据汇总.tsv            (panel b Han Gd-IgA1 P vs IgAN meta P)

写：
  阶段4/主图/out/Fig5_v3.1.{png,pdf,tiff}（投稿版；2 面板上下布局）

版本链：
  v1.0 2026-09-06 步骤 4.3 冻结（1×2 左右布局）
  v2.0 2026-09-06 重叠修复（v1.0 件留档）
  v3.0 2026-09-07 用户过审轮（_fig5_check.py + _diag5.py 双轨驱动）：
      (1) 布局重构：1×2 左右 → 2×1 上下（与 Fig1-4 全稿一致）。根因：
          左右布局时 panel a ax 仅 ~340 px 宽、x 尺度 ~100 px/unit，
          组内柱心距 0.30 unit = 30 px < 柱顶数字宽 45 px → EAS/EUR/meta
          相邻数字必叠 15 px OVERLAP。上下布局后 panel a 独占全宽 ~750 px、
          x 尺度 ~308 px/unit，柱心距 54 px > 数字宽，全部放下；
      (2) panel a 红框降矮：高 1.10 → 0.44（-0.04~0.40）让出顶部空间给
          legend；REVERSED 标签 va=top 贴框底下方；
      (3) panel a legend：fig-level (0.42, 0.955) 错位压标题 → 改 ax 内
          upper-right 图例 ncol=3 + frameon 白底；
      (4) panel b 标签错开：rs7856182 (x=8.62, y=1.00) 与 rs10238682
          (x=8.92, y=0.65) 同高文字 bbox 水平叠 72 px → rs7856182 标签改
          dy=0.75 放点上方，与 rs10238682 (dy=0.10) 纵向完全分离；
      (5) panel b y=x 标签 (8.5, 8.15) → (3.0, 3.7) 移到对角线上方空白；
      (6) panel b P=5e-8 标签 (0.5, 7.55) → (9.85, 7.45) 让位给 b 面板 legend；
      (7) panel b legend：fig-level (0.97, 0.475) 悬两面板间 → 改 ax 内
          upper-left ncol=1 + frameon 白底（水平虚线 y=7.3 穿过 → 白底遮）；
      (8) main() 移除两 fig.legend（与 ax 内 legend 重复且易错位）。
      自检 _fig5_check.py 五件套 TOTAL=0；_diag5.py 残留 <6px NEAR 全为
      同基线相邻柱顶数字/轴角假阳性（已目检确认无视觉粘连）。
  v3.1 2026-09-25 y 轴标签措辞对齐（配套标黄稿 v1.3）：
      panel a y 轴 "Effect-allele frequency (AF)" → "Allele frequency (AF)"。
      根因：正文 Fig. 5 图注已改为 "Allele frequencies"，Table 6 列名亦由
      "Effect allele" 改为 "Allele"（该表给的是"所列等位"的频率，仅 IgAN meta
      柱为 EAF），故图像标签须与图注/表头同口径。仅措辞变化、无几何改动
      （新标签更短，不引入重叠），panel a 图例 "IgAN meta EAF (GCST90018866)"
      保留（该柱确为 EAF）。旧件 Fig5_v3.0.* 留档不删。
      （注：panel_b 内联 "v3.1：" 注释系 v3.0 冻结前的内部迭代标签，其改动
      均已包含在 v3.0 冻结件中，与本次 v3.1 发布无关。）
  v3.7 2026-10-08 用户目检修正（本次）：panel b legend 下移，让出水平虚线。
      根因：legend 白底（frameon=True, framealpha=0.95）以 y=0.995 框位
      压在 -log10(5e-8)=7.30 的水平虚线上，实测该虚线行橙色像素在
      x 274→758 px（约 484 px）整段被遮。改 bbox_to_anchor y 0.995→0.85，
      legend 顶边由数据 y≈8.33 降到 ≈7.10，落在虚线下方，虚线全线贯通。
      仅几何改动，无数据/数值/结论变化。
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
    OKABE_ITO, C_NOTTEST, C_GW_LINE,
    FS_PANEL_LABEL, FS_AXIS_LABEL, FS_TICK, FS_LEGEND, FS_INCELL, FS_TITLE,
    FIG_WIDTH_MM, DPI_DRAFT, save_final, mm_to_in,
)

PROJECT_ROOT = _paths.LEGACY
F5A_FREQ = PROJECT_ROOT("阶段3/M5_东亚/M5_gnomAD_EASEUR_频率_20260903.tsv")
F5B_EVD  = PROJECT_ROOT("阶段3/M5_东亚/lead_x_东亚证据汇总.tsv")
OUT_DIR  = PROJECT_ROOT("阶段4/主图/out")
os.makedirs(OUT_DIR, exist_ok=True)


def panel_a(ax):
    """a. 4 lead × 3 人群等位基因频率对比（EAS/NFE 柱为 AF；IgAN meta 柱为 EAF）"""
    df = pd.read_csv(F5A_FREQ, sep="\t")
    leads = ["rs13226913", "rs10238682", "rs7856182", "rs5910940"]

    x = np.arange(len(leads))
    # v3.0：组内柱心距 0.25 → 0.30、柱宽 0.25 → 0.20（柱间 0.10 间隙）。
    # 根因：原 0.25/0.25 零间隙紧贴，rs5910940 组 EAS(0.546)/EUR(0.511)
    # 柱高几乎相等 → 两柱顶数字在同一水平线，中心距 0.25 unit < 文字
    # 半宽和 → bbox 重叠 24×1px（自检 TXT）。扩心距后 0.30 unit ≈ 30 px
    # > 5 字符数字半宽和 ~22 px，稳。
    step = 0.30
    width = 0.20

    eas = [float(df[df["lead_rsid"] == r].iloc[0]["gnomAD_EAS_AF"]) for r in leads]
    eur = [float(df[df["lead_rsid"] == r].iloc[0]["gnomAD_EUR_NFE_AF"]) for r in leads]
    meta = []
    for r in leads:
        v = df[df["lead_rsid"] == r].iloc[0]["IgAN_meta_EAF"]
        meta.append(float(v) if pd.notna(v) else np.nan)

    # 三柱：EAS (红) + EUR (蓝) + meta (灰)
    ax.bar(x - step, eas, width=width,
           color=OKABE_ITO["vermillion"], edgecolor="black", lw=0.4,
           label="EAS (gnomAD v4)", zorder=3)
    ax.bar(x, eur, width=width,
           color=OKABE_ITO["blue"], edgecolor="black", lw=0.4,
           label="EUR/NFE (gnomAD v4)", zorder=3)
    for i, m in enumerate(meta):
        if not np.isnan(m):
            ax.bar(x[i] + step, m, width=width,
                   color=OKABE_ITO["gray"], edgecolor="black", lw=0.4,
                   label="IgAN meta EAF (GCST90018866)" if i == 0 else None,
                   zorder=3)
            ax.text(x[i] + step, m + 0.012, f"{m:.3f}", color="black",
                    fontsize=FS_TICK - 1, ha="center", va="bottom")
        else:
            # rs5910940 X chr NA：灰格
            ax.add_patch(plt.Rectangle((x[i] + step - width/2, 0), width, 0.04,
                                       facecolor="#EEEEEE", edgecolor=C_NOTTEST,
                                       hatch="///", lw=0.5, zorder=2))
            ax.text(x[i] + step, 0.08, "chrX\nNA*",
                    color=C_NOTTEST, fontsize=FS_TICK - 1,
                    ha="center", va="bottom", style="italic")

    for i, (e, u) in enumerate(zip(eas, eur)):
        ax.text(x[i] - step, e + 0.012, f"{e:.3f}", color=OKABE_ITO["vermillion"],
                fontsize=FS_TICK - 1, ha="center", va="bottom")
        ax.text(x[i], u + 0.012, f"{u:.3f}", color=OKABE_ITO["blue"],
                fontsize=FS_TICK - 1, ha="center", va="bottom")

    # v3.6 2026-10-08：按审稿意见删除 rs7856182 的"反差红框 + REVERSED"标注。
    # 根因：把该位点标为"REVERSED / 方向相反"与图注、正文的新口径
    #（"is the least common in East Asians"——强调"东亚频率最低"，
    # 而非"方向相反"）不一致；红框还会被读成"该位点行为异常"。
    # 现仅保留三柱本身（EAS / EUR / IgAN meta），不加任何高亮或"reversed"语义。
    # 连带：下方 ylim 下界由 -0.28 收到 -0.06（原为给 REVERSED 标签留位）。

    ax.set_xticks(x)
    # X tick 两行（rs ID 上、gene 下）+ y 标签让位（v3.6：REVERSED 标签已删，
    # 下界原为 -0.28 的留位随之收紧到 -0.06）。
    ax.set_xticklabels([f"$\\it{{{l}}}$\n{g}" for l, g in zip(
        leads, ["C1GALT1", "C1GALT1", "GALNT12", "C1GALT1C1"])],
        fontsize=FS_TICK - 0.5)
    ax.tick_params(axis="x", pad=2)
    ax.set_ylabel("Allele frequency (AF)", fontsize=FS_AXIS_LABEL)
    ax.set_ylim(-0.06, 1.20)
    ax.tick_params(axis="y", labelsize=FS_TICK)
    # ax-level 图例改为 fig-level（不被 constrained_layout 强制回内）
    panel_a.legend_handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=OKABE_ITO["vermillion"], edgecolor="black", lw=0.4),
        plt.Rectangle((0, 0), 1, 1, facecolor=OKABE_ITO["blue"], edgecolor="black", lw=0.4),
        plt.Rectangle((0, 0), 1, 1, facecolor=OKABE_ITO["gray"], edgecolor="black", lw=0.4),
    ]
    panel_a.legend_labels = ["EAS (gnomAD v4)", "EUR/NFE (gnomAD v4)", "IgAN meta EAF (GCST90018866)"]
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    ax.set_title("a  Lead × pop. AF",
                 fontsize=FS_PANEL_LABEL, loc="left", pad=4, weight="bold")
    # v3.0：图例放 ax 右上图内（红框已降矮腾出顶部 y>0.96 空白区）。
    # 三柱语义色（EAS 红/EUR 蓝/meta 灰）+ ncol=3 单行。放 ax 内而非
    # fig-level：panel a 独占全宽后右上空白足够。
    ax.legend(handles=panel_a.legend_handles, labels=panel_a.legend_labels,
              loc="upper right", bbox_to_anchor=(0.995, 0.995),
              ncol=3, frameon=True, framealpha=0.95, edgecolor="#CCCCCC",
              fontsize=FS_LEGEND - 0.5, borderaxespad=0.3)


def panel_b(ax):
    """b. Han Gd-IgA1 -log10P vs IgAN meta -log10P 散点"""
    df = pd.read_csv(F5B_EVD, sep="\t")
    leads = ["rs13226913", "rs10238682", "rs7856182", "rs5910940"]

    x_vals = []
    y_vals = []
    labels = []
    colors = []
    sizes = []
    rids = []
    for r in leads:
        row = df[df["lead_rsid"] == r].iloc[0]
        pg = float(row["Wang2021_Han_GdIgA1_P"])
        pi = row["igan_meta_P"]
        x_vals.append(-np.log10(pg))
        rids.append(r)
        if pd.isna(pi) or str(pi).strip() in ("", "nan"):
            y_vals.append(None)
            labels.append(f"{r}\n(X chr NA)")
            colors.append(OKABE_ITO["bluish_grn"])
            sizes.append(140)
        else:
            pi_f = float(pi)
            y_vals.append(-np.log10(pi_f))
            labels.append(r)
            colors.append(OKABE_ITO["vermillion"] if pg < 1e-4 else OKABE_ITO["blue"])
            sizes.append(170 if pg < 1e-4 else 100)

    # 标签偏移：rs10238682 (x=8.92, y=0.65) 与 rs7856182 (x=8.62, y=1.00)
    # 两点相邻，若都放 (x+0.2, y+0.1) va=bottom 则两文字 bbox 水平叠
    # 72px（_diag5.py 实测）。rs7856182 标签改放点上方 dy=0.75，与
    # rs10238682 标签 (y≈0.75~0.93) 纵向完全错开（1.75 起）。
    for x, y, r, lbl, c, sz in zip(x_vals, y_vals, rids, labels, colors, sizes):
        dy = 0.75 if r == "rs7856182" else 0.10
        if y is not None:
            ax.scatter(x, y, s=sz, c=c, edgecolors="black", lw=1.0, zorder=4)
            ax.text(x + 0.2, y + dy, lbl, color="black",
                    fontsize=FS_TICK - 0.5, ha="left", va="bottom", weight="bold")
        else:
            ax.scatter(x, 1.5, s=sz, c=c, edgecolors="black", lw=1.0,
                       marker="^", zorder=4)
            ax.text(x + 0.2, 1.5, lbl, color=c,
                    fontsize=FS_TICK - 0.5, ha="left", va="center", weight="bold")

    lim_max = 10.0
    ax.plot([0, lim_max], [0, lim_max], color="gray", lw=0.8, ls=":", zorder=2)
    # v3.0：y=x 标签从 (lim_max*0.95, lim_max*0.95+0.3)=(9.5, 9.8) 移入
    # ylim 0-8.5 内 (8.5, 8.15)——原位置 y=9.8 超出 ylim 8.5，文字画在
    # ax 顶缘之外（clip 关闭时溢出压标题区/被 tight 裁切）。
    # v3.1：再移到 (3.0, 3.7)——对角线 y=x 在 x=3.0~3.5 处 y≤3.5 < 3.7，
    # 文字完全在对角线之上方空白区（x 2.5~3.6 无任何散点/标签）。
    ax.text(3.0, 3.7, "y = x", color="gray",
            fontsize=FS_TICK - 0.5, ha="left", va="bottom", style="italic")

    ax.axvline(-np.log10(5e-8), color=C_GW_LINE, lw=1.0, ls="--", zorder=2)
    ax.axhline(-np.log10(5e-8), color=C_GW_LINE, lw=1.0, ls="--", zorder=2)
    # v3.1：P=5e-8 标签从左上 (0.5, 7.55) 移到右上 (9.85, 7.45)——左上
    # 需让位给 b 面板 legend（upper-left）；右上 x 8.9~9.9 无散点
    # （红/蓝点 y≤1.5、散点标签 y≤1.95），水平虚线 y=7.3 下方贴线。
    ax.text(9.85, 7.45, "P=5e-8", color=C_GW_LINE,
            fontsize=FS_TICK - 0.5, ha="right", va="bottom")

    ax.set_xlabel(r"Wang'21 Han Gd-IgA1 $-\log_{10}(P)$",
                  fontsize=FS_AXIS_LABEL)
    ax.set_ylabel(r"IgAN meta (Sakaue'21) $-\log_{10}(P)$",
                  fontsize=FS_AXIS_LABEL)
    ax.set_xlim(0, lim_max)
    ax.set_ylim(0, 8.5)
    ax.tick_params(axis="both", labelsize=FS_TICK)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    # v3.1：b 面板 legend 从 fig-level 移入 ax 内左上（数据区 x0~2.5 无点，
    # 唯 P=5e-8 标签已右移腾位；水平虚线 y=7.3 穿过 → frameon=True 白底
    # 遮线保证可读）。ncol=1 三行：P<1e-4 红 / nominal 蓝 / X chr 绿三角。
    # v3.7（2026-10-08，用户目检）：legend 顶边原在 y=0.995 框位 = 数据
    # y≈8.33，白底正好压住水平虚线 y=7.30（-log10(5e-8)）在 x=0.05~2.46
    # 一段（实测该行橙色像素 263→2273 px 中间断口 274→758 px，共 ~484 px
    # 虚线被遮）。改为 bbox_to_anchor y=0.85 → legend 顶边降到数据 y≈7.10，
    # 整个框落在虚线下方，虚线全线贯通；框底 ≈ y 5.37，仍远高于该区最近的
    # 数据点（rs13226913 蓝点 y=0.5）与 "y = x" 标签（x=3.0 在框右缘之外），
    # 不引入新重叠。
    ax.legend(handles=[plt.scatter([], [], s=170, c=OKABE_ITO["vermillion"],
                                   edgecolors="black", lw=1.0,
                                   label="Han Gd-IgA1 P<1e-4"),
                       plt.scatter([], [], s=100, c=OKABE_ITO["blue"],
                                   edgecolors="black", lw=1.0,
                                   label="Han Gd-IgA1 nominal"),
                       plt.scatter([], [], s=140, c=OKABE_ITO["bluish_grn"],
                                   edgecolors="black", lw=1.0, marker="^",
                                   label="X chr (no IgAN)")],
              loc="upper left", bbox_to_anchor=(0.005, 0.85),
              ncol=1, frameon=True, framealpha=0.95, edgecolor="#CCCCCC",
              fontsize=FS_LEGEND - 0.5, borderaxespad=0.4)
    ax.set_title("b  Gd-IgA1 vs IgAN P",
                 fontsize=FS_PANEL_LABEL, loc="left", pad=4, weight="bold")


def main():
    # v3.0：1×2 左右布局 → 2×1 上下布局（与 Fig1-4 全稿一致 + 根治
    # panel a 数字重叠）。根因：左右布局时 panel a ax 仅 ~340 px 宽、
    # x 尺度 ~100 px/unit，组内柱心距 0.30 unit=30 px < 柱顶数字宽
    # 45 px → EAS/EUR/meta 相邻数字必叠（_diag5.py 实测 15px 宽
    # OVERLAP + 多组 NEAR<6px）。a 独占全行后 ax ~750 px、x 尺度
    # ~180 px/unit，柱心距 0.30 unit=54 px > 45 px，全部数字放下。
    fig = plt.figure(figsize=(mm_to_in(220), mm_to_in(205)), dpi=DPI_DRAFT)
    from matplotlib.gridspec import GridSpec
    gs = GridSpec(2, 1, figure=fig, height_ratios=[0.92, 1.0], hspace=0.42)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[1, 0])

    panel_a(ax_a)
    panel_b(ax_b)

    # v3.1：legend 全部改为 panel_a/panel_b 函数内 ax.legend（a 右上图内 /
    # b 左下图内）。此前 fig-level 双 legend（(0.42,0.955) 压 panel a 标题、
    # (0.97,0.475) 悬在两面板空隙）与 ax.legend 重复且易错位，全部移除。

    out_stem = os.path.join(OUT_DIR, os.environ.get("OUT_STEM", "Fig5_v3.1"))
    save_final(fig, out_stem)
    plt.close(fig)
    print("WROTE", out_stem + ".png")


if __name__ == "__main__":
    main()