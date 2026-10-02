"""fig1_render.py — Fig.1 概念框架（2026-09-06 步骤 4.3 制作第七图）
作者：AI 王严课题组代办｜依赖：matplotlib/pandas/numpy（fig43 venv）

读：
  阶段1/M1_GEO表达锚点/GSE73953_IgAN_vs_HC_best.csv       (panel b PBMC)
  阶段1/M1_GEO表达锚点/GSE115857_IgAN_vs_Ctrl_LD_best.csv  (panel b 肾 bulk)
  阶段1/M1_GEO表达锚点/GSE93798_IgAN_vs_Ctrl_best.csv      (panel b 肾小球)

写：
  阶段4/主图/out/Fig1_v4.3.{png,pdf,tiff}（投稿版；3 垂直面板）

版本链：
  v1.0  2026-09-06 步骤 4.3 初版（3 垂直面板：a 概念 + b 热图 + c 证据链）
  v1.1  2026-09-06 整合审计升版（v1.0→v1.1 仅文案/字号微调）
  v2.0  2026-09-06 重叠修复（panel a/c 字号/列宽/布局整体调整）
  v3.0  2026-09-07 反馈再修复（panel a second hit 注释下移避开 Gd-IgA1 节点框；
        "no shared causal structure" 标签下移避开底部双层文字；
        panel c 列宽重分配 + Theme 文字精简 + Key conclusion 字号下调至
        7.5pt，确保 Theme/Key conclusion/Q/Fig 四列文字完全不重叠）。
  v3.1  2026-09-07 反馈再修复 v2（panel a second hit 框高 1.2→1.4 让三行文字
        完全在框内、no shared 标签同步上移 0.1 unit；panel c Key conclusion
        列宽 4.90→4.55+文字精简 47~60 字符到 45~55 字符，确保单行不溢出
        压住 Q/Fig 列；Theme 列宽 2.50→2.45 配合 Key 起点前移）。
  v3.2  2026-09-07 Key conclusion 字号 7.5pt→7pt（FS_TICK-1.5），实测
        7.5pt 字符宽 ≈ 0.065 inch/字符，50 字符仍溢出 4.55 列宽压住 Q 列；
        7pt 字符宽 ≈ 0.048 inch/字符，55 字符单行占 ~3.4 unit，4.55 unit
        列宽留 1.1 unit 余量，彻底消除与 Q/Fig 列的紧贴/重叠。
  v3.3  2026-09-07 反馈再修复 v3（面板 a/c 大改：节点 title 贴框顶、body 居中
        下部，title/body 间隔 ≥0.8 unit；表头蓝带与首行彻底分离；7 行 L0–L6
        完整可见。自检新增 → TOTAL 0 text/patch overlap）。
  v3.4  2026-09-07 反馈再修复 v4：扩大所有框以彻底容纳文字（用户反馈"字
        母都到框框外面了"）。
        • 节点框 2.6×2.9 → **3.1×2.95**；x 位置 1.5/5.0/8.5 → 1.7/5.0/8.3；
          节点 3 title 改双行 "Mesangial IgA / → IgAN"；body 节点 2 缩短为
          "Galactose-deficient IgA1 / (Wang'21 P=1.2e-9)"；body 节点 3 改为
          3 行 "(risk axis null) / rs13226913 P=0.28 / rs10238682 P=0.22"。
        • second hit 框 2.2×1.2 → **3.1×1.65**（高 +0.45 让 3 行 8pt 文字
          完全在框内）。
        • 底部双层 4.5×2.0 → **4.5×2.10**（高 +0.10）；矩形下移 0.1 给 body
          留底 padding。
        • panel c 表头蓝带 0.45 → **0.60 unit**；header_y 10.575 → 11.30；
          ylim 扩 0-11 → 0-12 给表头留 0.40 unit 顶 padding；GridSpec
          height_ratios [100,62,96] → [100,62,108]；fig 高度 272→280mm。
        • 自检脚本 _fig1_check.py 新增 check_containment：文字中心是否在
          所属 patch 内、patch 内是否完全包住文字 bbox（tol 3px），自动
          列出所有 ESCAPE；v3.4 自检必须 TOTAL=0 ESCAPE 才通过。
  v3.7  2026-09-10 图文一致性修复（阶段4 ↔ 阶段4.5 全量对齐）：
        • panel c L4 行 Fig 列 "Tbl.S3" → "Tbl.S4"——L4 直接查证实际
          交付于补充表 **S4a/S4b**（S3 是 SuSiE 信号级复检表），原指代
          为陈旧编号；
        • panel c L6 行 "5-azaC [20]" → 全拼 "5-azacytidine [20]"——
          "5-azaC" 是正文、图注与附录缩写表均未使用的孤儿缩写，同图
          panel a 与 Fig.6 均已全拼，属同图内不一致；
        • panel a/c 的 max PPH4 统一为三位小数 **0.031**（原 a="0.03"、
          c="0.0308"），与 Fig.1 图注 "maximum PPH4 = 0.031" 及 Results
          正文 "maximum 0.031" 对齐；
        • panel c L7 行 Fig 列 "Fig.S5" → "Fig.S5-6"——L7 上游调控层的结论
          同时交付于补充图 **S5**（scRNA 细胞图谱）与 **S6**（CollecTRI
          regulon 活性），Fig. 1 图注 (c) 亦写 "Figures S5–S6 and Tables
          S9–S10"。
        未改动任何数据、版式几何或结论。
  v3.8  2026-09-11 重渲染件已交付（Fig1_v3.8.*），**该版未在本版本链留条目**
        （docstring 停在 v3.7）；产物留档于 `阶段4/主图/out/`，可随时重渲。
        ——本轮据实补记，不回填来源不明的描述。
  v3.9  2026-09-25 panel c L6 行引文编号对齐 x3 docx 文献表（配套标黄稿 v1.5）：
        "C1GALT1C1=6.0 > C1GALT1=4.5; DNMTi **[19, 20]**" →
        "…DNMTi **[22, 23]**"。
        根因：图件原按**母本线**（29 条文献）渲染，该线中 [19]=Lin 2018、
        [20]=Sun 2015（即 DNMTi 证据：IL-17→5-azacytidine 恢复表达；
        Cosmc 启动子甲基化→地西他滨逆转）。而投稿用的 **x3 docx 线**
        （33 条）里 [19]=Wu 2013、[20]=Dong 2021 是 **C1GALT1 促瘤**文献，
        DNMTi 证据已移到 **[22]=Lin 2018 / [23]=Sun 2015**；x3 稿的 Fig. 6
        图注亦已写 [22]/[23]。**同一张图图文编号体系不一致**，故改画面编号
        与之对齐（Fig. 6 panel d 同步：`fig6_render.py` 升 v3.6）。
        仅一处字符串变化，无几何/数据改动；新串比旧串长 3 字符，L6 行
        "Key conclusion" 列实测仍不越界（`_fig1_check.py` 全 0）。
        旧件 Fig1_v3.8.* 留档不删。
  v4.0  2026-09-25 panel c L6 行研究问题编号对齐正文（图文一致性修复，配套标黄稿 v1.7）：
        L6（Intervention）Q 列 "Q4" → **空**。
        根因：x3 稿 Background 将 **Q4 定义为“效应是否具人群特异性”**，
        A-7 轮已删除挂在干预小节（Results 第 7 节）的 "(Q4)"，即 Q4 只属于
        L5（EAS heterogeneity）。本图 L5/L6 都标 Q4，沿用早期
        `fig1_caption_v01.md` 的旧定义“干预窗与人群特异性”→ L6 的 Q4 与正文不符。
        干预层是全文四问的**落点/综合**（Background："so as to identify the regulatory
        layer at which … interventions would be most tractable"），不对应任何单一 Q，故清空该格。
        仅改一个单元格字符串，col_w 与行几何不变。旧件 Fig1_v3.9.* 留档不删。
  v4.1  2026-09-27 九层结构转正（配套标黄稿 v1.13）：
        第 9 个 Results 小节（单细胞细胞状态跨检查；Methods 已命名
        "Single-cell B-cell reanalysis (independent cell-level cross-check)"）
        由"八层链之外的独立跨检查"正式转为**第九层 L8**——新增 L8 行
        （Theme "Cell-state check"；Key conclusion 取实测锚点
        "C1GALT1C1 detection 3.1→5.2% (naive→plasma)"（Results R8：检出率
        3.1%→5.2%，Fisher P = 2.8e-25）；Q 列 Q1；Fig 列 Fig.S4，对应
        Additional file 1: Figure S4 与 Tables S7/S8）。
        面板标题 "(8 layer conclusions)" → "(9 layer conclusions)"、L0–L7 → L0–L8。
        为在 ylim 0–12 内排下 9 行，行距 1.115 → 1.10（框高 0.94 不变）：
        末行 L8 框底 y = 0.63；表头底 10.925 与首行框顶 10.37 的 0.555 unit
        间距不变；相邻框间空隙 0.175 → 0.16 unit。仅新增一行，无数据、无
        坐标重排、无结论改动。旧件 Fig1_v4.0.* 留档不删。
  v4.2  2026-09-27 行序与 Results 小节序统一（配套标黄稿 v1.14）：
        第 8/9 行**对调**——原第 8 行 "L7 Upstream regulon"、第 9 行
        "L8 Cell-state check"，与正文次序相反（Results 第 8 小节＝细胞状态
        跨检查，第 9 小节＝上游调控层；第 9 小节以 "Within the same atlas"
        起句，依赖第 8 小节引入的 GSE285335 图谱，故正文次序不可动）。
        统一方向＝**图随正文**：新第 8 行＝("L7","Cell-state check"，
        "C1GALT1C1 detection 3.1→5.2% (naive→plasma)"，Q1，Fig.S4)、
        新第 9 行＝("L8","Upstream regulon"，"Cytokine-DNMT1-Sp/KLF
        circuit; emp. P<0.05"，Q2，Fig.S5-6)。
        连带：面板标题 L0–L8 不变（层数仍为 9），ch8 图注 "the eighth/
        ninth row" 两句对调、ch1 摘要 Methods 枚举末两项对调。
        **仅行位置（含各行 Q/Fig 单元随之移动）改变，无文字、无数据、无
        几何参数改动**（行距 1.10、框高 0.94、列宽、字号均不变；行底色
        奇偶填充随位置重排，属版式副作用）。旧件 Fig1_v4.1.* 留档不删。
  v4.3  2026-10-01 上游调控层（L8）改阴性结论（方案 B；配套 ch4 v2.6 / ch8 v2.7）：
        该层改按**完整 CollecTRI 网络**（64,516 条记录；表达面板 4,199→8,248 基因；
        778 regulons）并按**各调控子自身靶数**匹配零分布重算 → 无任一调控子过
        FDR 后仍超其匹配零分布，该层**整体为阴性**。panel c 第 9 行 Key conclusion
        "Cytokine-DNMT1-Sp/KLF circuit; emp. P<0.05" →
        **"Negative vs matched null; Sp/KLF-dense island"**（该层仅存的
        **序列层**观察：C1GALT1C1 启动子岛密集分布可被甲基化阻断的 Sp/KLF 基序）。
        仅此一处字符串变化：第 9 行仍然存在、层数仍为 9、面板标题 "L0–L8 evidence
        chain (9 layer conclusions)" 不变、第 8 行（L7 Cell-state check）不动、
        行距/框高/列宽/字号/坐标一律不变，无数据改动。新串比旧串短 2 字符，
        L8 行 "Key conclusion" 列实测仍不越界（`_fig1_check.py` 全 0）。
        旧件 Fig1_v4.2.* 留档不删。
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
from matplotlib.patches import FancyBboxPatch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from style_common import (  # noqa: E402
    OKABE_ITO, C_NOTTEST, C_GW_LINE,
    FS_PANEL_LABEL, FS_AXIS_LABEL, FS_TICK, FS_LEGEND, FS_INCELL, FS_TITLE,
    FIG_WIDTH_MM, DPI_DRAFT, save_final, mm_to_in,
)

PROJECT_ROOT = _paths.LEGACY
F1B_73953 = PROJECT_ROOT("阶段1/M1_GEO表达锚点/GSE73953_IgAN_vs_HC_best.csv")
F1B_115857 = PROJECT_ROOT("阶段1/M1_GEO表达锚点/GSE115857_IgAN_vs_Ctrl_LD_best.csv")
F1B_93798 = PROJECT_ROOT("阶段1/M1_GEO表达锚点/GSE93798_IgAN_vs_Ctrl_best.csv")
OUT_DIR  = PROJECT_ROOT("阶段4/主图/out")
os.makedirs(OUT_DIR, exist_ok=True)


def new_fig_three_panels_v(panel_h_mm=(80, 70, 95), w_mm=FIG_WIDTH_MM):
    """3 垂直面板（专为 Fig1 设计）。"""
    total_h_mm = sum(panel_h_mm) + 14
    fig, axes = plt.subplots(
        3, 1,
        figsize=(mm_to_in(w_mm), mm_to_in(total_h_mm)),
        dpi=DPI_DRAFT,
        constrained_layout=True,
        gridspec_kw={"height_ratios": panel_h_mm, "hspace": 0.45},
    )
    return fig, axes


def panel_a(ax):
    """a. 概念框架（horizontal flow）— v3.4 反馈再修复：扩大所有框以
    完全容纳文字。

    v3.3 自检仍检出文字溢出框外（patch#1/#2/#3/#4/#5 共 10 处 ESCAPE）：
      patch#1 Gd-IgA1 body "Galactose-deficient...Wang'21 1.2e-9"
             长 ~30 字符 8pt → ~1.86 inch，框 2.6 unit ≈ 1.64 inch 不足
      patch#2 Mesangial title "Mesangial IgA → IgAN" 19 字符 11pt bold
             → ~1.95 inch，框 1.64 inch 不足；body 38 字符 → 1.9 inch
      patch#3 second hit 3 行 8pt 字符总高 ~ 0.45 inch，框 1.2 unit
             ≈ 0.42 inch 不足
      patch#4/#5 底部 body "DOES NOT DRIVE IgAN RISK (coloc 14/25...)"
             "ACTIONABLE WINDOW (C1GALT1C1...)" 2 行 8pt，框 2.0 unit
             ≈ 0.7 inch 临界

    v3.4 修复（按用户要求"扩大框"为主，文字适度精简为辅）：
      节点框：2.6×2.9 → **3.1×2.95**（宽 +0.5，高 +0.05）；
        x 位置 1.5/5.0/8.5 → 1.7/5.0/8.3（让节点 3 不再顶右边缘）；
        节点 3 title 改双行 "Mesangial IgA / → IgAN" 防溢出；
        body 节点 2 缩短为 "Galactose-deficient IgA1 / (Wang'21 P=1.2e-9)"；
        body 节点 3 改为 3 行 "(risk axis null) / rs13226913 P=0.28 /
        rs10238682 P=0.22"。
      second hit 框：2.2×1.2 → **3.1×1.65**（宽 +0.9，高 +0.45）；
        位置 (3.45, 4.5)，文字 3 行居中。
      底部双层：4.5×2.0 → **4.5×2.10**（高 +0.10）；
        矩形位置下移 0.1（0.5→0.4），title/body 不变以保持视觉一致。
    """
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 11)
    ax.axis("off")

    # ---- 顶部三节点（v3.5：节点 3 title 单行 10pt + body 2 行解决 v3.4 框内
    # title/body 重叠；其他节点保持 v3.4 双行 title + 短 body） ----
    nodes = [
        (1.7, 8.3, "B / plasma cells",
         ["IgA1 O-glycosylation axis", "(5 genes)"],
         "#E5F1F7", FS_TITLE),
        (5.0, 8.3, "Gd-IgA1 ↑",
         ["Galactose-deficient IgA1", "(Wang'21 P=1.2e-9)"],
         "#FFE9C9", FS_TITLE),
        # 节点 3：title 单行（19 字符 11pt bold ≈ 1.99 inch 临界）→ 字号降到
        # FS_TITLE-1=10pt（19 chars × 0.092 ≈ 1.75 inch < 框宽 1.95）；
        # body 2 行（不再 3 行），避免与 title 互相重叠
        (8.3, 8.3, "Mesangial IgA → IgAN",
         ["(risk null: rs13226913 P=0.28,", "rs10238682 P=0.22)"],
         "#FFD6D6", FS_TITLE - 1),
    ]
    for x, y, title, body_lines, bg, title_fs in nodes:
        # 框宽 3.1, 高 2.95（容纳 title 多行 + body 多行 + padding）
        box = FancyBboxPatch((x - 1.55, y - 1.475), 3.1, 2.95,
                              boxstyle="round,pad=0.05",
                              facecolor=bg, edgecolor="black", lw=1.5)
        ax.add_patch(box)
        # title：贴框顶 y+0.85（距框顶 y+1.475 留 0.6 unit padding）
        ax.text(x, y + 0.85, title, color="black",
                fontsize=title_fs, ha="center", va="center", weight="bold")
        # body：合并多行用 \n；居中偏下
        body_text = "\n".join(body_lines)
        body_cy = y - 0.60
        ax.text(x, body_cy, body_text, color="black",
                fontsize=FS_TICK - 0.3, ha="center", va="center")

    # 节点间实箭头（gap 仅 0.20 units ≈ 0.13 inch ~ 39px@300dpi，箭头仍清晰）
    arrow_kw = dict(arrowstyle="->", color="black", lw=1.5, mutation_scale=16)
    for x1, x2 in [(3.25, 3.45), (6.55, 6.75)]:
        ax.annotate("", xy=(x2, 8.3), xytext=(x1, 8.3), arrowprops=arrow_kw)

    # ---- "second hit" 注释（v3.5：框 3.4×1.65，"mucosal immunity / complement"
    # 31 chars 8pt ≈ 1.92 inch ≈ 3.05 unit，框 3.4 unit 余 0.35 unit padding） ----
    box = FancyBboxPatch((3.30, 4.5), 3.4, 1.65,
                          boxstyle="round,pad=0.05",
                          facecolor="white", edgecolor=C_GW_LINE, lw=1.0)
    ax.add_patch(box)
    # 2 行 8pt 字符总高 ~ 0.26 inch ~ 0.74 unit；框 1.65 unit 上下各留 0.45 unit padding
    ax.text(5.0, 5.32, "second hit:\nmucosal immunity / complement",
            color=C_GW_LINE, fontsize=FS_TICK - 0.5, ha="center", va="center",
            style="italic")
    # 虚线箭头 second hit → Mesangial 框底（绕过 Gd-IgA1 框下方空区）
    ax.annotate("", xy=(8.3, 6.85), xytext=(5.0, 6.15),
                arrowprops=dict(arrowstyle="->", color=C_GW_LINE, lw=1.0,
                                linestyle="--", mutation_scale=14,
                                connectionstyle="arc3,rad=-0.20"))

    # ---- 底部双层（v3.4：框 4.5×2.10，下移 0.1 给 body 留底 padding） ----
    # Germline genetics（左）
    ax.add_patch(plt.Rectangle((0.3, 0.4), 4.5, 2.10,
                               facecolor="#FFEEEE", edgecolor="black", lw=1.0))
    ax.text(2.55, 2.05, "Germline genetics",
            color="#D55E00", fontsize=FS_TITLE, ha="center", va="center", weight="bold")
    ax.text(2.55, 1.10,
            "DOES NOT DRIVE IgAN RISK\n(coloc 14/25 null; max PPH4 = 0.031)",
            color="black", fontsize=FS_TICK - 0.5, ha="center", va="center")
    # Acquired epigenetic（右）
    ax.add_patch(plt.Rectangle((5.2, 0.4), 4.5, 2.10,
                               facecolor="#E5F5E5", edgecolor="black", lw=1.0))
    ax.text(7.45, 2.05, "Acquired epigenetic",
            color="#00755E", fontsize=FS_TITLE, ha="center", va="center", weight="bold")
    ax.text(7.45, 1.10,
            "ACTIONABLE WINDOW\n(C1GALT1C1 promoter; 5-aza DNMTi)",
            color="black", fontsize=FS_TICK - 0.5, ha="center", va="center")

    # ---- "no shared causal structure"标签：位于底部双层顶(2.5)与 second hit 底(4.5)空带 ----
    ax.text(5.0, 3.5, "no shared\ncausal structure",
            color="black", fontsize=FS_TICK - 0.3, ha="center", va="center",
            style="italic",
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                      edgecolor="black", lw=0.5))

    ax.set_title("a  Conceptual framework of the axis",
                 fontsize=FS_PANEL_LABEL, loc="left", pad=4, weight="bold")


def panel_b(ax):
    """b. 5 基因 × 3 数据集 logFC 方向性"""
    genes = ["C1GALT1", "C1GALT1C1", "GALNT12", "GALNT2", "ST6GALNAC2"]

    df_pbmc = pd.read_csv(F1B_73953, sep=",").set_index("gene")
    df_bulk = pd.read_csv(F1B_115857, sep=",").set_index("gene")
    df_glomer = pd.read_csv(F1B_93798, sep=",").set_index("gene")

    logfcs = np.zeros((3, len(genes)))
    pvals = np.zeros((3, len(genes)))
    for i, g in enumerate(genes):
        logfcs[0, i] = float(df_pbmc.loc[g, "logFC"])
        pvals[0, i] = float(df_pbmc.loc[g, "pvalue"])
        logfcs[1, i] = float(df_bulk.loc[g, "logFC"])
        pvals[1, i] = float(df_bulk.loc[g, "pvalue"])
        logfcs[2, i] = float(df_glomer.loc[g, "logFC"])
        pvals[2, i] = float(df_glomer.loc[g, "pvalue"])

    cmap_colors = ["#2874A6", "#AED6F1", "#FFFFFF", "#F5B7B1", "#A93226"]
    from matplotlib.colors import ListedColormap, BoundaryNorm
    cmap = ListedColormap(cmap_colors)
    norm = BoundaryNorm([-5, -1.5, -0.2, 0.2, 1.5, 5], cmap.N)

    sig_mask = (pvals < 0.05) & (np.abs(logfcs) > 1.0)

    im = ax.imshow(logfcs, cmap=cmap, norm=norm, aspect="auto")

    for i in range(3):
        for j in range(len(genes)):
            txt = f"{logfcs[i, j]:+.2f}"
            if sig_mask[i, j]:
                txt += "*"
            color = "white" if abs(logfcs[i, j]) > 1.5 else "black"
            ax.text(j, i, txt, ha="center", va="center",
                    color=color, fontsize=FS_TICK - 0.5)

    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels([f"$\\it{{{g}}}$" for g in genes],
                       fontsize=FS_TICK, rotation=0)
    ax.set_yticks(range(3))
    # Y 标签精简：缩写 + n（与正文数据集口径一致，避免左边界裁切）
    ax.set_yticklabels(
        ["GSE73953 (PBMC)\nn = 15/2 pooled",
         "GSE115857 (kidney bulk)\nn = 55/7, DASL",
         "GSE93798 (glomerular)\nn = 20/22"],
        fontsize=FS_TICK - 1)
    ax.tick_params(axis="y", pad=2)
    ax.set_title("b  Expression direction across 3 datasets (logFC)",
                 fontsize=FS_PANEL_LABEL, loc="left", pad=4, weight="bold")
    cbar = ax.figure.colorbar(im, ax=ax, fraction=0.04, pad=0.04, shrink=0.85)
    cbar.ax.tick_params(labelsize=FS_TICK - 0.5)
    cbar.set_label("logFC (IgAN vs Ctrl)", fontsize=FS_TICK - 0.5)

    # 星号注释：放 axes 下方独立 fig-level 区域（不与 colorbar 冲突）
    ax.text(0.0, -0.32,
            "* P<0.05 & |logFC|>1; PBMC uniformly down (5/5); kidney directionally separated",
            transform=ax.transAxes, fontsize=FS_TICK - 0.5,
            ha="left", va="top", color="black")


def panel_c(ax):
    """c. L0–L8 证据链条带（表格）— v4.3（2026-10-01）：第 9 行 L8 改阴性结论
      （"Negative vs matched null; Sp/KLF-dense island"），层数/几何不变。
      v4.2（2026-09-27）：第 8/9 行对调，
      使行序与 Results 小节序一致（第 8 行＝细胞状态跨检查 L7、第 9 行＝
      上游调控层 L8）；行数仍为 9、几何不变。
      v4.1（2026-09-27）：行数 8→9，新增 L8
      “Cell-state check（单细胞细胞状态跨检查）”行（原为八层链之外的独立跨检查，
      v1.13 起正式为第九层）；行距 1.115→1.10 以在 ylim 0–12 内排下 9 行，
      框高 0.94 不变。
      v3.6（2026-09-10）：行数 7→8，新增 L7
      “Upstream regulon（cytokine–DNMT1–Sp/KLF circuit）”行以配合证据层 L7
      整合；行距 1.32→1.115、框高 1.00→0.94 以在 ylim 0–12 内排下 8 行。
      v3.4 反馈再修复（保留）：
      (1) v3.3 表头蓝带高 0.45 unit 实测仍被 9pt 表头字号顶/底各顶出 4~5px
          （patch#0~4 ESCAPE 共 12 处）→ 蓝带高 0.45→0.60 unit，header_y 上移
          至 11.30（ylim 扩到 0-12 给表头留 0.70 unit padding）；
      (2) v3.3 列宽/数据行/字号不变（已自检通过）。
    列宽：Layer 0.15~0.80 / Theme 0.85~3.05 / Key 3.10~8.55 / Q 8.60~9.10 /
    Fig 9.10~9.85。
    """
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 12)
    ax.axis("off")

    rows = [
        ("L0", "Expression",
         "5 genes down in PBMC (GSE73953; -1.7 to -4.4)",
         "Q1", "Fig1b"),
        ("L1", "Cell-type cis-eQTL",
         "C1GALT1 B_naive P=9.4e-6; GTEx P=5.6e-12",
         "Q1", "Fig2"),
        ("L2", "Coloc with IgAN",
         "Coloc 14/25 OK; max PPH4=0.031; alleles null",
         "Q3", "Fig3"),
        ("L3", "Methylation",
         "5 gene-body CpG (no promoter); MR null",
         "Q2", "Fig4"),
        ("L4", "cis-eQTL lookup",
         "10/10 axis-region top variants not IgAN-sig.",
         "Q3", "Tbl.S4"),
        ("L5", "EAS heterogeneity",
         "rs13226913 EAS 0.93 vs EUR 0.42; rs7856182 inv.",
         "Q4", "Fig5"),
        ("L6", "Intervention",
         "C1GALT1C1=6.0 > C1GALT1=4.5; DNMTi [22, 23]",
         "", "Fig6"),
        # v4.2（2026-09-27）：第 8/9 行对调 → 行序与 Results 小节序一致
        # （正文第 8 小节＝细胞状态跨检查、第 9 小节＝上游调控层）
        ("L7", "Cell-state check",
         "C1GALT1C1 detection 3.1→5.2% (naive→plasma)",
         "Q1", "Fig.S4"),
        # v4.3（2026-10-01）：L8 上游调控层改为**阴性**——完整 CollecTRI 网络 +
        # 按各调控子自身靶数匹配零分布重算后无任一调控子过 FDR 仍超其匹配零分布；
        # 该层仅存的观察是**序列层**的 C1GALT1C1 启动子岛 Sp/KLF 基序密度。
        ("L8", "Upstream regulon",
         "Negative vs matched null; Sp/KLF-dense island",
         "Q2", "Fig.S5-6"),
    ]

    headers = ["Layer", "Theme", "Key conclusion", "Q", "Fig"]
    # v3.7 列宽复核（图文一致性修复）：
    #   • Theme 列 2.15 → 2.40 unit：v3.7 实测 "EAS heterogeneity"（17 字符
    #     8.7pt bold ≈ 2.15 unit）正好顶满原列宽，与 Key 列文字视觉连读
    #     （渲染图上读作 "EAS heterogeneityrs13226913"）；加宽后 Theme/Key
    #     之间余 0.25 unit 空白。
    #   • Figs 列（第 5 列）0.90 → 0.65 unit：文字中心由 9.50 左移至 9.375，
    #     使最宽的 "Fig.S5-6"（8 字符 8pt 斜体 ≈ 0.50 unit）右缘 9.76 落在
    #     行框右缘 9.85 之内，`_fig1_check.py` ESCAPE=0。
    #   • Key 列 5.10 → 4.80 unit：仍够容纳最长 47 chars × 0.048 in = 2.26 in
    #     ≈ 3.59 unit 的 Key 文字 + 右侧余量。
    col_x = [0.10, 1.05, 3.50, 8.30, 9.05]
    col_w = [0.95, 2.40, 4.80, 0.75, 0.65]

    # 表头蓝带：v3.5 高 0.75 unit（容纳 9pt 表头字符 + 上下各 0.12 unit padding）；
    # header_y=11.30，带范围 10.925~11.675（ax 顶 12 → 表头顶 11.675 余 0.325 unit）
    header_y = 11.30
    for h, x, w in zip(headers, col_x, col_w):
        ax.add_patch(plt.Rectangle((x, header_y - 0.375), w, 0.75,
                                   facecolor="#A8C8DC", lw=0.8))
        ax.text(x + w / 2, header_y, h, color="black",
                fontsize=FS_TICK, ha="center", va="center", weight="bold")

    # v4.1（2026-09-27）：加入 L8 细胞状态跨检查层 → 9 数据行。y0=9.90，行距 1.10，
    # 框高 0.94 → 末行 L8 y=1.10（框底 0.63，仍在 ylim 内）；表头底 10.925 与首行
    # 框顶 10.37 间隔 0.555 unit（未变）；相邻框间空隙 0.16 unit。
    # v3.6（2026-09-10）：加入 L7 上游调控层 → 8 数据行。y0=9.90，行距 1.115，
    # 框高 0.94 → 末行 L7 y=2.095（框底 1.625，距 ax 底余量充足）；表头底
    # 10.925 与首行框顶 10.37 间隔 0.555 unit；相邻框间空隙 0.175 unit（≈4.2pt）。
    # 出图后以 _fig1_check.py 复核无 text/patch 越界与重叠。
    for i, (layer, theme, conclusion, q, fig) in enumerate(rows):
        y = 9.90 - i * 1.10
        bg = "#F5F5F5" if i % 2 == 0 else "white"
        ax.add_patch(plt.Rectangle((0.15, y - 0.47), 9.70, 0.94,
                                   facecolor=bg, edgecolor="gray", lw=0.4))
        ax.text(col_x[0] + col_w[0] / 2, y, layer, color=C_GW_LINE,
                fontsize=FS_TICK, ha="center", va="center", weight="bold")
        ax.text(col_x[1] + 0.05, y, theme, color="black",
                fontsize=FS_TICK - 0.3, ha="left", va="center", weight="bold")
        # Key conclusion：7pt 单行 ≤46 字符（约 4.0 unit 宽），列宽 5.45 充足
        ax.text(col_x[2] + 0.05, y, conclusion, color="black",
                fontsize=FS_TICK - 1.5, ha="left", va="center")
        ax.text(col_x[3] + col_w[3] / 2, y, q, color="black",
                fontsize=FS_TICK - 1.0, ha="center", va="center")
        ax.text(col_x[4] + col_w[4] / 2, y, fig, color="black",
                fontsize=FS_TICK - 1.0, ha="center", va="center", style="italic")

    ax.set_title("c  L0–L8 evidence chain (9 layer conclusions)",
                 fontsize=FS_PANEL_LABEL, loc="left", pad=4, weight="bold")


def main():
    fig = plt.figure(figsize=(mm_to_in(183), mm_to_in(280)), dpi=DPI_DRAFT)
    from matplotlib.gridspec import GridSpec
    gs = GridSpec(3, 1, figure=fig, height_ratios=[100, 62, 108],
                  hspace=0.42)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[2, 0])

    panel_a(ax_a)
    panel_b(ax_b)
    panel_c(ax_c)

    out_stem = os.path.join(OUT_DIR, os.environ.get("OUT_STEM", "Fig1_v4.3"))
    save_final(fig, out_stem)
    plt.close(fig)
    print("WROTE", out_stem + ".png")


if __name__ == "__main__":
    main()