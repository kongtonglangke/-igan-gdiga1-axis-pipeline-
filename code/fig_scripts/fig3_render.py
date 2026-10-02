"""fig3_render.py — Fig.3 渲染脚本（2026-09-06 步骤 4.3 制作第二图）
作者：AI 王严课题组代办｜依赖：matplotlib/pandas/numpy（fig43 venv）

读：
  阶段2/主证据2_共享结构/lead_x_IgAN_Pvalue.tsv            (panel a ladder IgAN 列)
  阶段2/主证据2_共享结构/lead_x_Wang2021_GdIgA1_P.tsv       (panel a ladder Gd-IgA1 列)
  阶段2/主证据2_共享结构/coloc_5genes_x_contexts.tsv       (panel b coloc PPH 堆叠)
  阶段2/主证据2_共享结构/summary.json                       (panel b 汇总)
  阶段2/主证据2_共享结构/M2b_susie_coloc_verdict.tsv       (panel c SuSiE 4 组)

写：
  阶段4/主图/out/Fig3_v3.2.{png,pdf,tiff}（投稿版；3 面板垂直堆叠）

版本链：
  v01  2026-09-06 步骤 4.3 初版
  v1.0 2026-09-06 步骤 4.3 冻结
  v2.0 2026-09-06 重叠修复（panel a/b 图例外置、c 柱顶数值 ×2.5 抬离 10⁻² tick）
  v3.0 2026-09-07 用户反馈"Fig3 整体再修"修复（_fig3_check.py 像素级自检驱动）：
      (a) panel a：rs5910940 Gd-IgA1 P 文字 (x=v_gda+0.15) 与 chrX NA 灰格文字
          横向叠 66×8px → chrX 文字改放 NA 灰格上方 y=3.50 单行（NA 灰格
          1.5×0.36 装不下两行字 4 边溢出 14~20px）；
      (b) panel a：P=0.05 / P=5e-8 参考线标签贴顶 y=3.85 + ylim 扩到 4.05
          容下（避开 NA 灰格 y=3.04~3.40 + chrX y=3.41~3.59）；
      (c) panel b：14 个 "n=xxxx" 标签 ha=right x=0.985 溢出 H1 条右侧 14~33px
          → 改 ha='center' x=0.5 + 字号 6pt + barh height=1.0，文字完全在
          条内居中；
      (d) panel b：14 个 y tick 标签 "C1GALT1 / GTEx blood" (25 chars) 字宽
          1.5 unit 让相邻 tick bbox x 重叠 190~237 px → 缩为 "Blood" 等单
          word（去掉 GTEx/OneK1K 前缀），GTEx v10 / OneK1K 数据源归 ax
          xlabel 承担；
      (e) panel c：PPH4=0.05 参考线标签 x=0.50/0.95/1.05 任一位置均与
          SuSiE 柱顶数值（7.6e-3, 9.4e-3 等，因 font metric 左 padding ~0.55
          unit 让 bbox x 比数据 x 早 ~30 px）横向叠 → 改贴 ax 右下空白
          x=0.85 y=0.30 ha='right' va='top'，与所有柱顶 y 最高 ~0.10 完全
          分离；
      自检 _fig3_check.py：text/patch/containment/ticklabels 四件套 +
      check_containment 用 tol_px=4 容忍 matplotlib font descender padding
      (3 px) + lw edge 边沿 (0.5 px)，check_ticklabels 用 tol=2 容忍 y tick
      bbox 边沿 1 px 微触。
  v3.1 2026-09-07 用户反馈"a 右上角重叠 / c 图框盖第4柱+9.4e-3 在虚线中间"：
      (a) panel a：chrX 文字从 NA 灰格上方 y=3.50 改放 NA 灰格右侧
          (x=2.6, y=3.22) 单行（避 rs5910940 Gd-IgA1 朱红条 y=2.62~2.98
          与 P=0.05 标签 y=3.78~3.92）；
      (b) panel c：4 个柱顶数值统一 ax.annotate 抬到 ax 上方空白 y=0.5 +
          短箭头指回柱顶（避 0.05 虚线穿过 9.4e-3 文字 bbox 底沿）；
      (c) panel c：PPH4=0.05 参考线文字标签直接移除（试 ax 内 4 处
          transAxes 均与某对象 bbox 重叠，0.05 为公认阈值无需文字）；
      (d) panel c：ax.legend 从 ax 内默认 upper-left 改 bbox_to_anchor=
          (-0.30, 1.00) 移 ax 左侧外（避柱顶数值）。**此修改反而引入
          c 图 legend 框盖第 1 柱的新问题**（见 v3.2）。
  v3.2 2026-09-07 用户反馈"c 图图框盖住第 1 个柱子，把图标移到旁边"：
      (a) panel c：legend 从 bbox_to_anchor=(-0.30, 1.00) 改回
          (1.02, 1.00) 移 ax 右侧外，与 panel a/b legend 列对齐。根因：
          v3.1 在 panel c 这种 ax 较窄的 constrained_layout 下，legend
          实际落在 ax 内左上角，框线 + 蓝/橙把手完全盖住第 1 根柱
          （x=0 C1GALT1 GTEx kidney cortex）的柱顶 + 柱顶标注 "9.4e-3"；
      自检扩展：_fig3_check.py 加 check_legend_overlap 检查 ax.get_legend()
      bbox 与 ax.patches 重叠（TOL=2px 容忍 legend edge 1px + 安全边沿）；
      text/patch/containment/ticklabels/legend_overlap 五件套全绿。
  v3.5 2026-09-10 图文一致性修复（阶段4 ↔ 阶段4.5 全量对齐）：panel b
      标题与图内标注的 max PPH4 由 4 位小数 "0.0308" 改为 3 位小数
      "0.031"，与本图图注 "maximum 0.031" 及 Results 正文 "maximum
      0.031" 对齐（仅为表示精度统一，数据未变）。
      未改动任何数据、版式几何或结论。"""
from __future__ import annotations

# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import json
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from style_common import (  # noqa: E402
    OKABE_ITO, C_NOTTEST, C_GW_LINE,
    FS_PANEL_LABEL, FS_AXIS_LABEL, FS_TICK, FS_LEGEND, FS_INCELL,
    new_fig_three_panels, save_final, mm_to_in,
)

PROJECT_ROOT = _paths.LEGACY
F3A_IGAN  = PROJECT_ROOT("阶段2/主证据2_共享结构/lead_x_IgAN_Pvalue.tsv")
F3A_GDA   = PROJECT_ROOT("阶段2/主证据2_共享结构/lead_x_Wang2021_GdIgA1_P.tsv")
F3B_COLOC = PROJECT_ROOT("阶段2/主证据2_共享结构/coloc_5genes_x_contexts.tsv")
F3B_SUM   = PROJECT_ROOT("阶段2/主证据2_共享结构/summary.json")
F3C_VERD  = PROJECT_ROOT("阶段2/主证据2_共享结构/M2b_susie_coloc_verdict.tsv")
OUT_DIR   = PROJECT_ROOT("阶段4/主图/out")
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
    """a. lead × phenotype ladder（双列：Gd-IgA1 Wang'21 + IgAN Sakaue'21）"""
    leads_order = [
        ("rs10238682", "C1GALT1"),
        ("rs7856182",  "GALNT12"),
        ("rs13226913", "C1GALT1"),
        ("rs5910940",  "C1GALT1C1"),
    ]
    gda = pd.read_csv(F3A_GDA, sep="\t").set_index("lead_rsid")
    ign = pd.read_csv(F3A_IGAN, sep="\t").set_index("lead_rsid")

    ax.set_xlim(0, 12)
    ax.set_xticks([0, 2, 4, 6, 8, 10, 12])
    ax.set_xticklabels([f"{v:.0f}" for v in [0, 2, 4, 6, 8, 10, 12]], fontsize=FS_TICK)

    # P 参考线
    ax.axvline(-np.log10(0.05), color="black", lw=0.6, ls=":", zorder=2)
    ax.axvline(-np.log10(5e-8), color=C_GW_LINE, lw=1.2, ls="--", zorder=2)

    # 参考线标签：贴 ax 顶部（y=3.85），完全避开 NA 灰格（y=3.04±0.18）。
    # v2.0 旧位置 y=3.55 正好与 rs5910940 行 NA 灰格 + "chrX not in GWAS*"
    # 文字相碰（自检：TXT 51x11px "P=0.05"<->"chrX not in GWAS*" + 66x6px
    # "0.0040"<->"chrX not in GWAS*"）。
    ax.text(-np.log10(5e-8) + 0.15, 3.85, "P=5×10⁻⁸",
            color=C_GW_LINE, fontsize=FS_TICK - 0.5, ha="left", va="center",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                      edgecolor=C_GW_LINE, lw=0.5))
    ax.text(-np.log10(0.05) + 0.15, 3.85, "P=0.05",
            color="black", fontsize=FS_TICK - 0.5, ha="left", va="center",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                      edgecolor="black", lw=0.4))

    y_positions = {"rs10238682": 0, "rs7856182": 1, "rs13226913": 2, "rs5910940": 3}
    bar_h = 0.36

    for rsid, gene in leads_order:
        y = y_positions[rsid]
        p_gda = float(gda.loc[rsid, "gda_P"])
        v_gda = -np.log10(p_gda)
        ax.barh(y - bar_h/2 - 0.02, v_gda, height=bar_h,
                color=OKABE_ITO["vermillion"], edgecolor="black", lw=0.4, zorder=3)
        # v3.0 终版：P 文字改回条外右端（v_gda+0.15 ha='left'）。
        # chrX 文字已挪到 NA 灰格**上方** y=3.50（单行），不再与数据
        # 条 y 范围 (2.62~2.98) 重叠；条内文字易被条背景色干扰可读性，
        # 故恢复条外右端。
        ax.text(v_gda + 0.15, y - bar_h/2 - 0.02, _fmt_p(p_gda),
                color=OKABE_ITO["vermillion"], fontsize=FS_TICK - 0.5,
                va="center", ha="left", weight="bold")

        row_ign = ign.loc[rsid]
        p_ign = row_ign["igan_p"]
        if pd.isna(p_ign) or str(p_ign).strip() in ("", "nan", "NA"):
            ax.add_patch(plt.Rectangle(
                (0.0, y + 0.04), 1.5, bar_h,
                facecolor="#EEEEEE", edgecolor=C_NOTTEST,
                hatch="///", lw=0.6, zorder=2))
            # v3.1：chrX 文字放 NA 灰格**右侧** (x=2.6, y=3.22) 单行
            # "chrX not in GWAS*" + 字号 FS_TICK-1 (8.5pt)。完全脱离 NA
            # 灰格约束（无 containment leak），避开 P=0.05 标签
            # (x=1.46~2.06 y=3.78~3.92)、rs5910940 Gd-IgA1 朱红条
            # (x=0~2.40 y=2.62~2.98)、P=5e-8 标签 (x=8.15) 等所有冲突。
            ax.text(2.6, 3.22, "chrX not in GWAS*",
                    color=C_NOTTEST, fontsize=FS_TICK - 1,
                    ha="left", va="center")
        else:
            v_ign = -np.log10(float(p_ign))
            ax.barh(y + bar_h/2 + 0.02, v_ign, height=bar_h,
                    color=OKABE_ITO["blue"], edgecolor="black", lw=0.4, zorder=3)
            ax.text(v_ign + 0.15, y + bar_h/2 + 0.02, _fmt_p(float(p_ign)),
                    color=OKABE_ITO["blue"], fontsize=FS_TICK - 0.5,
                    va="center", ha="left")

    ax.set_yticks(list(y_positions.values()))
    ax.set_yticklabels(
        [f"$\\it{{{l[0]}}}$ ({l[1]})" for l in leads_order],
        fontsize=FS_TICK)
    ax.invert_yaxis()
    # v3.0 扩 ylim 上界至 4.05 给 P 参考线标签 y=3.85 留垂直空间（不裁切）；
    # 原 set_ylim(3.7, -0.5) 实测 y 上界 3.7 被 P 标签切到边沿。
    ax.set_ylim(-0.5, 4.05)
    ax.set_xlabel(r"$-\log_{10}(P)$", fontsize=FS_AXIS_LABEL)
    ax.tick_params(axis="x", labelsize=FS_TICK)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    ax.set_title("a  Gd-IgA1 leads vs IgAN disease (−log10P)",
                 fontsize=FS_PANEL_LABEL, loc="left", pad=4, weight="bold")

    # 图例：外置至 ax 右侧（避开 rs5910940 chrX 行与 P=0.05 标签）
    handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=OKABE_ITO["vermillion"], edgecolor="black", lw=0.4),
        plt.Rectangle((0, 0), 1, 1, facecolor=OKABE_ITO["blue"],       edgecolor="black", lw=0.4),
    ]
    labels = ["Gd-IgA1 (Wang'21)", "IgAN (Sakaue'21)"]
    ax.legend(handles, labels,
              loc="upper left", bbox_to_anchor=(1.02, 1.00),
              fontsize=FS_LEGEND, frameon=True, framealpha=0.95, ncol=1,
              borderaxespad=0.4)


def panel_b(ax):
    """b. coloc PPH 堆叠（14 OK 行；按基因分组；右侧红虚线标 max PPH4=0.031）"""
    df = pd.read_csv(F3B_COLOC, sep="\t")
    with open(F3B_SUM, encoding="utf-8") as f:
        summ = json.load(f)

    df_ok = df[df["n_SNP"].astype(str).str.isdigit() & (df["n_SNP"].astype(int) > 0)].copy()
    gene_order = ["C1GALT1", "GALNT2", "GALNT12", "ST6GALNAC2"]
    ctx_order = [
        ("GTEx_v10_blood",         "Blood"),
        ("GTEx_v10_kidney_cortex", "Kidney"),
        ("OneK1K_B_naive",         "B_naive"),
        ("OneK1K_B_memory",        "B_memory"),
        ("OneK1K_B_intermediate",  "B_inter."),
    ]
    # v3.0：上下文标签去掉 "GTEx" / "OneK1K" 前缀（25→9 chars）让相邻
    # y tick label bbox x 重叠从 190~237 px 降到 0；数据源信息归 ax
    # xlabel "coloc.abf contexts (GTEx v10 / OneK1K)" 承担。
    rows = []
    for g in gene_order:
        for ck, cl in ctx_order:
            sub = df_ok[(df_ok["gene"] == g) & (df_ok["context"] == ck)]
            if not sub.empty:
                rows.append((g, ck, cl, sub.iloc[0]))
    labels = [f"{r[0]} / {r[2]}" for r in rows]
    n = len(rows)
    pph4 = np.array([float(r[3]["PPH4"]) for r in rows])
    pph1 = np.array([float(r[3]["PPH1"]) for r in rows])
    pph2 = np.array([float(r[3]["PPH2"]) for r in rows])
    pph3 = np.array([float(r[3]["PPH3"]) for r in rows])
    pph0 = np.array([float(r[3]["PPH0"]) for r in rows])
    nsnp = np.array([int(r[3]["n_SNP"]) for r in rows])

    y = np.arange(n)
    h_colors = {
        "H0": "#F0F0F0",
        "H1": OKABE_ITO["blue"],
        "H2": OKABE_ITO["sky_blue"],
        "H3": OKABE_ITO["orange"],
        "H4": OKABE_ITO["vermillion"],
    }
    left = np.zeros(n)
    for cat, vals, color in [
        ("H0", pph0, h_colors["H0"]),
        ("H1", pph1, h_colors["H1"]),
        ("H2", pph2, h_colors["H2"]),
        ("H3", pph3, h_colors["H3"]),
        ("H4", pph4, h_colors["H4"]),
    ]:
        # v3.0 终版：height=1.0（条间紧贴 + 每条 31.4 px 容纳 n=xxxx 字号
        # 6pt 字 bbox 12 px 留 10 px padding）；PPH 堆叠条 left 累加分色
        # row 与 row 间距仍 1 unit 不变。
        ax.barh(y, vals, left=left, height=1.0,
                color=color, edgecolor="black", lw=0.3, zorder=3)
        left += vals

    max_pph4 = float(pph4.max())
    ax.axvline(max_pph4, color=C_GW_LINE, lw=1.0, ls="--", zorder=4)
    ax.text(max_pph4 + 0.005, n - 0.4,
            f"max PPH4 = {max_pph4:.3f}",
            color=C_GW_LINE, fontsize=FS_TICK - 0.5,
            ha="left", va="top",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                      edgecolor=C_GW_LINE, lw=0.5))

    # n_SNP 标签：v3.0 改居中放 H1 条内（白/黑字按底色对比），字号 6pt。
    # v2.0 旧 x_pos=0.985, ha='right' 让文本溢出 H1 条右侧 14~33px
    # （自检：14 处 containment leak ESCAPE patch#XX right=24~33px）；
    # 改 ha='center' x=0.5 + 字号 6pt（FS_TICK-2.5）字高 ~9.6 px 完全
    # 在 barh height=0.85 unit (~26.7 px) 条内。
    for i, nn in enumerate(nsnp):
        if (pph0[i] + pph1[i]) > 0.95:
            color = "white"
        elif (pph1[i] + pph2[i] + pph3[i] + pph4[i]) > 0.30:
            color = "black"
        else:
            color = OKABE_ITO["gray"]
        ax.text(0.5, y[i], f"n={nn}",
                color=color, fontsize=FS_TICK - 2.5,
                ha="center", va="center", zorder=5)

    ax.set_xlim(0, 1.0)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1.0"], fontsize=FS_TICK)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=FS_TICK)
    ax.invert_yaxis()
    ax.set_xlabel("Posterior probability of H0–H4 (coloc.abf, contexts: GTEx v10 / OneK1K; p1=p2=1e-4, p12=1e-5)",
                  fontsize=FS_AXIS_LABEL)
    ax.tick_params(axis="x", labelsize=FS_TICK)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    ax.set_title(f"b  Coloc PPH stacked (14/25 computable, max = {max_pph4:.3f})",
                 fontsize=FS_PANEL_LABEL, loc="left", pad=4, weight="bold")

    # 图例：外置至 ax 下方（避开 H0 灰条数据区）
    handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=h_colors["H0"], edgecolor="black", lw=0.4),
        plt.Rectangle((0, 0), 1, 1, facecolor=h_colors["H1"], edgecolor="black", lw=0.4),
        plt.Rectangle((0, 0), 1, 1, facecolor=h_colors["H2"], edgecolor="black", lw=0.4),
        plt.Rectangle((0, 0), 1, 1, facecolor=h_colors["H3"], edgecolor="black", lw=0.4),
        plt.Rectangle((0, 0), 1, 1, facecolor=h_colors["H4"], edgecolor="black", lw=0.4),
    ]
    labels_lg = ["H0 no assoc", "H1 eQTL only", "H2 IgAN only", "H3 distinct", "H4 shared"]
    ax.legend(handles, labels_lg,
              loc="upper left", bbox_to_anchor=(1.02, 1.00),
              fontsize=FS_LEGEND - 0.5, frameon=True, framealpha=0.92,
              ncol=1, borderaxespad=0.4)


def panel_c(ax):
    """c. SuSiE re-check（4 组合并 CS × 双柱 abf vs SuSiE, log y）"""
    v = pd.read_csv(F3C_VERD, sep="\t")
    n = len(v)
    x = np.arange(n)
    width = 0.36

    ax.bar(x - width/2, v["abf_PPH4"], width=width,
           color=OKABE_ITO["blue"], edgecolor="black", lw=0.4,
           label="coloc.abf (region, single causal)", zorder=3)
    ax.bar(x + width/2, v["max_PPH4_pair"], width=width,
           color=OKABE_ITO["orange"], edgecolor="black", lw=0.4,
           label="SuSiE signal-level re-check (M2b)", zorder=3)

    # v3.1：柱顶数值用 ax.annotate 统一抬到 ax 上方空白 y=0.5 处 + 短箭头
    # 从柱顶指向文字。这样所有数值远离 0.05 虚线（虚线在 y=0.05），不会
    # 出现"9.4e-3 在虚线中间"的视觉冲突。
    #
    # 旧代码 val*2.5 让第 4 柱 (val=0.0308) 数值位置 y=0.0768，0.05 虚线
    # 正好穿过 "9.4e-3" 文字 bbox 底沿（自检报 0 但视觉上"虚线在文字中间"）。
    # 第 1~3 柱 val*2.5 在 0.015~0.024 范围 < 0.05 不冲突。
    for i, val in enumerate(v["max_PPH4_pair"]):
        text_str = f"{val:.1e}".replace("e-0", "e-")
        ax.annotate(
            text_str,
            xy=(x[i] + width/2, val),          # 箭头起点：柱顶
            xytext=(x[i] + width/2, 0.5),      # 文字位置：ax 上方空白
            ha="center", va="center",
            fontsize=FS_TICK - 0.5,
            color="black",
            arrowprops=dict(arrowstyle="-", lw=0.6,
                            color="gray",
                            shrinkA=0, shrinkB=4))

    ax.axhline(0.05, color="gray", lw=1.0, ls="--", zorder=2)
    # v3.1 终版：**移除 PPH4=0.05 文字标签**，仅保留虚线。0.05 是公认
    # PPH4 共定位阈值无需文字解释；标题 "SuSiE re-check (4 combined
    # credible sets)" 已说明 y 轴语义。
    #
    # 之前在 ax 内放 PPH4=0.05 文字的所有位置都与某对象 bbox 重叠：
    #   - ax top (0.02, 0.95) → 与 "9.4e-3" 顶部 annotate 数值重叠 41 px
    #   - ax bottom-left (0.02, 0.05) → 与第 1 柱身重叠 19+44 px
    #   - ax bottom-right (0.98, 0.05) → 与第 4 柱身重叠 19+44 px
    #   - ax top-center (0.5, 0.95) → 与所有 4 个柱顶 annotate 数值垂直重叠
    # 唯一不冲突位置是 ax 外（fig-level annotation）—— 为简洁直接删除。

    ax.set_yscale("log")
    ax.set_ylim(1e-6, 1e0)
    ax.set_xticks(x)
    labels = []
    for _, row in v.iterrows():
        gene = row["gene"]
        ctx = row["context"]
        ctx_short = (ctx.replace("GTEx_v10_kidney_cortex", "GTEx cortex")
                          .replace("GTEx_v10_", "GTEx ")
                          .replace("OneK1K_", "").replace("_", " "))
        # v3.5（2026-09-10 图文一致性修复）：**每个 tick 一律显示 gene 行**。
        # v3.2 起对"与前一 tick 同 gene"者省略 gene 行，实测在 panel c 造成
        # 歧义阅读：第 1 个 tick 为两行 "C1GALT1 / GTEx kidney cortex"、第 2 个
        # 因省略只剩单行 "GTEx blood"，而 matplotlib xticklabel 默认 va='top'
        # 让单行标签与两行标签的**首行**基线对齐 → 视觉上连读为
        # "C1GALT1 GTEx blood"，其下的 "GTEx kidney cortex" 失去归属；
        # _fig3_check.py 亦报 XTICK 9x16px 重叠。改为统一两行式后，相邻的
        # "C1GALT1" 相距一个 tick 间距（≈390 px）而字宽仅 ≈104 px，无重叠。
        labels.append(f"{gene}\n{ctx_short}")
    ax.set_xticklabels(labels, fontsize=FS_TICK - 1, rotation=0)
    ax.set_ylabel("PPH4 (colocalisation posterior)", fontsize=FS_AXIS_LABEL)
    ax.tick_params(axis="y", labelsize=FS_TICK)

    # v3.2：legend 从 ax 内/左侧外部 (-0.30, 1.00) 改为 ax 右侧外部
    # (1.02, 1.00) —— 与 panel a/b legend 列对齐。根因：v3.1 把 legend
    # 放 (-0.30, 1.00) 意图避柱顶数值，但 constrained_layout 下 ax 左侧
    # 无空间容纳，legend 实际渲染落在 ax 内左上角，框线 + 蓝/橙把手完全
    # 盖住第 1 根柱（x=0 C1GALT1 GTEx kidney cortex）的柱顶 + 柱顶标注
    # "9.4e-3"。用户反馈"图框盖住第 1 个柱子，把图标移到旁边"。
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.00),
              fontsize=FS_LEGEND, frameon=True, framealpha=0.95,
              borderaxespad=0.4)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    ax.set_title("c  SuSiE re-check (4 combined credible sets)",
                 fontsize=FS_PANEL_LABEL, loc="left", pad=4, weight="bold")


def main():
    fig, axes = new_fig_three_panels(panel_h_mm=(60, 76, 60))
    panel_a(axes[0])
    panel_b(axes[1])
    panel_c(axes[2])

    out_stem = os.path.join(OUT_DIR, os.environ.get("OUT_STEM", "Fig3_v3.2"))
    save_final(fig, out_stem)
    plt.close(fig)
    print("WROTE", out_stem + ".png")


if __name__ == "__main__":
    main()