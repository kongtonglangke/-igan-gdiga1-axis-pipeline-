"""fig6_render.py — Fig.6 干预锚点（2026-09-06 步骤 4.3 制作第六图）
作者：AI 王严课题组代办｜依赖：matplotlib/pandas/numpy（fig43 venv）

读：
  阶段3/M6_干预/M6_可干预性评分.tsv            (panel a,b 评分；4 维 G/D/F/P)
  阶段3/M6_干预/M6_靶点药物证据表.tsv          (panel b 干预证据；panel c,d 用药细节)

写：
  阶段4/主图/out/Fig6_v3.6.{png,pdf,tiff}（投稿版；2×2 网格 a/b/c/d）

设计要点（v1.1 整合审计升版：GALNT2 配体数 3 research ligands、红色不可成药区
axvspan 覆盖 Tdark+Tbio、干预列 5-azacytidine、panel d 出口措辞对齐 R7）：
  - a 5 基因可干预性总分堆叠条（G/D/F/P）
  - b 基因 × 干预手段热图（DNMTi(5-azacytidine)/EGCG/HDACi class I/HDACi research/sialyltransferase）
  - c Pharos 直接成药性 × 总分（提示干预杠杆须经转录/表观层）
  - d Pharos→表观层出口示意（5-azacytidine 用词=正文 R7）

版本链：
  v1.0 2026-09-06 步骤 4.3 冻结
  v1.1 2026-09-06 整合审计升版（GALNT2 3 ligands、axvspan 覆盖、d 措辞对齐 R7）
  v2.0 2026-09-06 重叠修复整体版（v1.0/v1.1 件留档）
  v3.0 2026-09-07 用户过审轮（_fig6_check.py + _diag6.py 双轨驱动，13 处归零）：
      (c) panel c：ST6GALNAC2 标签 ha=left 右伸 166px 与 GALNT2 长标签头
          在 y 202~208 带水平叠 22px → 改 ha='right' 放点左侧 (dx=-0.08)；
      (d) panel d 整体重排（原 11pt 单行标题超框右溢 168px + title/body
          垂直仅 8px 叠）：改**分行堆叠引擎**——title_lines/body_lines
          列表 + 行 pitch（title 0.86 / body 0.72 unit）+ 内容块垂直居中 +
          动态 box 高（红 3.08 / 绿 3.94 / 深绿 3.94 unit）+ 从顶向下排布
          + ylim 0-14 + 箭头改 x=5 间隙内；title 8.5pt / body 7.2pt；
          ✓✓ 标记中心 x 8.45 防右溢 9.5；
      (b) panel b：x tick 标签 FS_TICK→FS_TICK-1（"HDACi res." 与相邻
          标签字面仅 1.5~1.7px，降 1pt 拉出 ~12px 间隙）。
      自检 _fig6_check.py 五件套 TOTAL=0（新增：跳过空文本 annotate 伪
      text、跳过 #FFCCCC 扣分斜纹配对、panel d set_xticks([]) 消不可见
      tick）；_diag6.py 邻近度残留 = panel c 标签横跨粉底带边界（Tbio/
      Tdark 语义本就在不可成药带内，设计意图）+ panel d 行间 3.4~3.7px
      （300dpi 成图 ~5px，视觉分离成立）。
  v3.1 2026-09-07 用户二轮反馈（"a 图的图标放到整张图最下面 / c 图
      圆点与字母重合 + 一串字母越界到 d 图"）：
      (a) panel a 图例位置：原 fig.legend bbox_to_anchor=(0.50, -0.02)
          把 G/D/F/-P 四色块孤悬整张图最下方中央，与 a 标题相距过远
          → 用户看不出"这是 a 图的图标"。改 ax 内右下
          (loc='lower right', bbox_to_anchor=(1.00, -0.05)) ncol=2
          frameon=False——色块紧贴 a 图 xlabel 上方，归属明确。
          main() 删除 fig.legend。
      (c) panel c 标签二轮修复：
          - 统一删 ligand 数（"3 research ligands"）：信息已在 x 轴副
            标签 "Tchem (3 ligands)" 表达，散点标签只保留 (Txxx) tier；
          - 5 点全按 y 错位避免水平叠（v3.0 同行 3 点 GALNT12/
            ST6GALNAC2/GALNT2 起点 x 互差 1 unit 但 ha='left' 全右伸
            导致标签 bbox 互相覆盖）：
            C1GALT1C1 (2, 6.0)  上右 dy=+0.55
            C1GALT1   (2, 4.5)  上右 dy=+0.55
            GALNT12   (1, 2.5)  上右 dy=+0.42
            ST6GALNAC2(2, 2.5)  下左 dy=-0.75（避 ytick "2" 行 1.95~2.05）
            GALNT2    (3, 2.5)  上左 dy=+0.95 ha='right'（避 panel c 右缘
            越界进入 d 图；与 GALNT12 垂直 0.53 unit 充分分离）。
      自检 _fig6_check.py 修正：check_containment 跳过 alpha<1 的
      axvspan patch（"Not directly druggable" 粉底默认占满 ax ylim
      高度，所有 text 中心 y 都在其内，bbox 越 axvspan 边界是设计意图
      不应报 ESCAPE）。五件套 TOTAL=0。
  v3.2–v3.5（2026-09-07 … 09-11）中间版已交付（Fig6_v3.2/3.3/3.4/3.5.*），
      **本 docstring 未逐版留条目**（链条停在 v3.1）；各版改动见包内 README
      与 `图表发表性审查报告_v1.0/v1.1`。本轮据实补记，不回填来源不明的描述。
  v3.6 2026-09-25 panel d 出口引文编号对齐 x3 docx 文献表（配套标黄稿 v1.5）：
      "reversible by decitabine **[20]**" → "reversible by decitabine **[23]**"。
      根因：图件原按**母本线**（29 条文献）渲染，该线 [20]=Sun 2015（Cosmc
      启动子甲基化→地西他滨逆转，即本行所指证据）；投稿用的 **x3 docx 线**
      （33 条）里 [20]=Dong 2021 是 **C1GALT1 促瘤**文献，该证据已移到
      **[23]=Sun 2015**，且 x3 稿的 Fig. 6 图注（(b) 段）亦已写 "[23]"。
      **同一张图图文编号体系不一致**，故改画面编号与之对齐（Fig. 1
      panel c L6 同步：`fig1_render.py` 升 v3.9）。
      仅一处字符串变化，无数据/几何改动；新串与旧串字符数相同（20→23 同为
      2 字符），panel d 行布局不变。旧件 Fig6_v3.5.* 留档不删。
  v3.7 2026-09-25 panel a x 轴标签对齐评分公式上限（图文一致性修复，配套标黄稿 v1.7）：
      "Intervenability score (0–10)" → "Intervenability score (0–8)"。
      根因：评分公式 G(0–3)+D(0–3)+F(0–2)−P(0–2) 的理论上限为 **8**；x3 稿
      Fig. 6 图注与 Table 7 表头（A-9 轮）已均改为 (0–8)，panel c 的 y 轴亦为
      0–8，**惟 panel a 的 x 轴标签仍写 0–10** → 同一张图内两套量程。
      仅改标签字符串，无数据/几何改动（`set_xlim(0, 9.5)` 保留为总分标注留白，
      matplotlib 自动刻度仍为 0/2/4/6/8，无 10 刻度）。旧件 Fig6_v3.6.* 留档不删。"""
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
F6_SCORE = PROJECT_ROOT("阶段3/M6_干预/M6_可干预性评分.tsv")
F6_TARG  = PROJECT_ROOT("阶段3/M6_干预/M6_靶点药物证据表.tsv")
OUT_DIR  = PROJECT_ROOT("阶段4/主图/out")
os.makedirs(OUT_DIR, exist_ok=True)


def panel_a(ax):
    """a. 5 基因 × 4 维堆叠条 (G/D/F/-P)"""
    df = pd.read_csv(F6_SCORE, sep="\t")
    # 排序：按总分降序
    df = df.sort_values("可干预性总分(0-8)", ascending=False).reset_index(drop=True)

    genes = df["基因"].tolist()
    G = df["评分_遗传锚点G(0-3)"].values
    D = df["评分_药物手柄D(0-3)"].values
    F = df["评分_功能证据F(0-2)"].values
    P = df["评分_风险扣分P(0-2)"].values
    total = df["可干预性总分(0-8)"].values

    x = np.arange(len(genes))
    # 堆叠（从左到右）：G → D → F → -P（红色斜纹为负向扣分）
    ax.barh(x, G, color=OKABE_ITO["blue"], edgecolor="black", lw=0.4,
            label="Genetic anchor G", zorder=3)
    ax.barh(x, D, left=G, color=OKABE_ITO["orange"], edgecolor="black", lw=0.4,
            label="Drug handle D", zorder=3)
    ax.barh(x, F, left=G+D, color=OKABE_ITO["bluish_grn"], edgecolor="black", lw=0.4,
            label="Functional evidence F", zorder=3)
    # 风险扣分用红色斜纹覆盖
    for i, p in enumerate(P):
        if p > 0:
            ax.barh(x[i], -p, left=G[i]+D[i]+F[i], height=0.6,
                    color="#FFCCCC", edgecolor=C_GW_LINE, lw=0.6,
                    hatch="///", label="Risk penalty -P" if i == 0 else None, zorder=4)
    # 总分文字
    for i, t in enumerate(total):
        ax.text(t + 0.15, i, f"{t:.1f}", color="black",
                fontsize=FS_TICK, ha="left", va="center", weight="bold")

    ax.set_yticks(x)
    # 基因名斜体（mathtext）
    italic_labels = ["$\\it{" + g + "}$" for g in genes]
    ax.set_yticklabels(italic_labels, fontsize=FS_TICK)
    ax.invert_yaxis()  # C1GALT1C1 最高分置顶（conventional）
    ax.set_xlim(0, 9.5)
    ax.set_xlabel("Intervenability score (0–8)", fontsize=FS_AXIS_LABEL)
    ax.tick_params(axis="both", labelsize=FS_TICK)
    # v3.2 修正：v3.1 注释声称改 ax.legend 但代码漏调用，渲染 a 图无图例。
    # 改用 fig-level legend（图级 legend 不受 constrained_layout 强制回
    # ax 内的限制），bbox_to_anchor 横向居中于 a 图正下方：
    # a 图在 GridSpec gs[0,0] 横向占 fig 左半 [0.05, 0.48]，纵向
    # [0.52, 0.96]（fig 180mm 高，a/c 共 1 行 = 90mm）。bbox_to_anchor=
    # (0.265, 0.495) 落在 a/c 间隙中央偏 a 底；4 色块 ncol=4 横向排
    # 列总宽约 60mm，frameon=False 视觉轻量。
    panel_a.legend_handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=OKABE_ITO["blue"], edgecolor="black", lw=0.4),
        plt.Rectangle((0, 0), 1, 1, facecolor=OKABE_ITO["orange"], edgecolor="black", lw=0.4),
        plt.Rectangle((0, 0), 1, 1, facecolor=OKABE_ITO["bluish_grn"], edgecolor="black", lw=0.4),
        plt.Rectangle((0, 0), 1, 1, facecolor="#FFCCCC", edgecolor=C_GW_LINE, lw=0.6, hatch="///"),
    ]
    panel_a.legend_labels = ["Genetic anchor G", "Drug handle D", "Functional evidence F", "Risk penalty −P"]
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    ax.set_title("a  Intervenability score (G/D/F/−P)", fontsize=FS_PANEL_LABEL,
                 loc="left", pad=4, weight="bold")


def panel_b(ax):
    """b. 基因 × 干预手段热图（证据等级 0/1/2/3）"""
    df = pd.read_csv(F6_TARG, sep="\t")
    # 干预手段列：5 类（精简 X 标签为单行 + 副说明括号）
    interventions = [
        ("DNMTi",        "5-aza class"),
        ("EGCG",         "non-nucleoside"),
        ("HDACi I",      "vorinostat"),
        ("HDACi res.",   "research"),
        ("Sialyl-Ti",    "inhibitor"),
    ]
    matrix = {
        "C1GALT1C1":   [3, 2, 1, 1, 0],
        "C1GALT1":     [2, 2, 2, 1, 0],
        "GALNT12":     [0, 1, 0, 0, 0],
        "GALNT2":      [0, 0, 0, 0, 0],
        "ST6GALNAC2":  [0, 0, 0, 0, 1],
    }
    genes = ["C1GALT1C1", "C1GALT1", "GALNT12", "GALNT2", "ST6GALNAC2"]
    mat = np.array([matrix[g] for g in genes])

    # 热图（基色：0=浅灰，1=浅蓝，2=蓝，3=深蓝）
    cmap_colors = ["#F5F5F5", "#D4E6F1", "#7FB3D5", "#2874A6"]
    from matplotlib.colors import ListedColormap
    cmap = ListedColormap(cmap_colors)
    im = ax.imshow(mat, cmap=cmap, vmin=-0.5, vmax=3.5, aspect="auto")
    # 注释数字
    for i in range(len(genes)):
        for j in range(len(interventions)):
            txt_color = "white" if mat[i, j] >= 2 else "black"
            ax.text(j, i, str(mat[i, j]), ha="center", va="center",
                    color=txt_color, fontsize=FS_TICK)

    ax.set_xticks(range(len(interventions)))
    # X 标签：精简为单行 + 副说明走 ax 下方 caption（避免水平叠加）。
    # v3.0：FS_TICK→FS_TICK-1——"HDACi res."(~10ch) 与相邻标签字面
    # 仅 1.5~1.7px（imshow 5 列 ~84px/列 vs 标签宽 ~82px），降 1pt 拉出
    # ~12px 间隙。
    # v3.3：FS_TICK-1→FS_TICK-2（7.5pt→6.5pt），"HDACi res." 9ch×6.5pt
    # ≈ 4.7mm < 5 列均宽 ~12.7mm + 单字画间距缩 ~30%，标签间横向
    # 视觉间距从 ~1.5px 增到 ~8px，彻底消除"EGCG 与 HDACi 紧贴"。
    ax.set_xticklabels([it[0] for it in interventions], fontsize=FS_TICK - 2)
    ax.tick_params(axis="x", pad=4)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels([f"$\\it{{{g}}}$" for g in genes], fontsize=FS_TICK)
    ax.set_title("b  Gene × intervention evidence (0–3)",
                 fontsize=FS_PANEL_LABEL, loc="left", pad=4, weight="bold")
    # 副说明（X 轴下方一行）+ 图例外置右上
    ax.set_xlabel("Interventions (5-aza class / EGCG / vorinostat / research / sialyl-Ti)",
                  fontsize=FS_TICK - 0.8, labelpad=4)
    legend_elements = [
        plt.Rectangle((0, 0), 1, 1, facecolor="#F5F5F5", edgecolor="black", lw=0.4,
                      label="0: none / not druggable"),
        plt.Rectangle((0, 0), 1, 1, facecolor="#D4E6F1", edgecolor="black", lw=0.4,
                      label="1: theoretically reachable"),
        plt.Rectangle((0, 0), 1, 1, facecolor="#7FB3D5", edgecolor="black", lw=0.4,
                      label="2: in vitro functional"),
        plt.Rectangle((0, 0), 1, 1, facecolor="#2874A6", edgecolor="black", lw=0.4,
                      label="3: direct IgAN/clinical"),
    ]
    ax.legend(handles=legend_elements, loc="upper left",
              bbox_to_anchor=(1.02, 1.00), fontsize=FS_LEGEND - 0.5,
              frameon=True, framealpha=0.95, ncol=1, borderaxespad=0.4)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)


def panel_c(ax):
    """c. Pharos 直接成药性 × 总分（散点）"""
    df = pd.read_csv(F6_SCORE, sep="\t")
    # Pharos tclass 映射（C1GALT1, C1GALT1C1, ST6GALNAC2 = Tbio; GALNT2 = Tchem; GALNT12 = Tdark）
    pharos = {
        "C1GALT1C1":  ("Tbio", 0),
        "C1GALT1":    ("Tbio", 0),
        "GALNT12":    ("Tdark", 0),
        "GALNT2":     ("Tchem", 3),
        "ST6GALNAC2": ("Tbio", 0),
    }
    tclass_x = {"Tdark": 1, "Tbio": 2, "Tchem": 3}

    # 标签布局：5 点按 y 错位避免水平叠。v3.2（用户反馈
    # "最右面那个圆没有字母了——让修改不是让删除"）：
    #   GALNT2 标签从 v3.1 (dx=-0.10, dy=+0.95) ha='right' 圆点
    #   左上方 1.0 unit 外——圆点旁看似无字。改回 (dx=+0.10, dy=+0.55)
    #   ha='left' 圆点正右上方 (3.10, 3.05)：实际文字 ~0.31 unit
    #   远在 xlim 3.7 内，水平 [3.10, 3.41] 与 GALNT12 [1.10, 1.41]
    #   间距 1.69+ unit 无水平叠；垂直 y=3.05 与 GALNT12 y=2.92 差
    #   0.13 unit（@257px/unit ≈ 33px）充分分离。
    label_cfg = {
        "C1GALT1C1":   (+0.12, +0.55, "left"),   # 上右
        "C1GALT1":     (+0.12, +0.55, "left"),   # 上右
        "GALNT12":     (+0.10, +0.42, "left"),   # 上右
        "ST6GALNAC2":  (-0.05, -0.75, "right"),  # 下左（避 ytick "2"）
        "GALNT2":      (+0.10, +0.55, "left"),   # 上右（v3.2 圆点正右上方）
    }

    for _, row in df.iterrows():
        g = row["基因"]
        tc, nl = pharos[g]
        score = row["可干预性总分(0-8)"]
        x = tclass_x[tc]
        ax.scatter(x, score, s=160, c=OKABE_ITO["vermillion"],
                   edgecolors="black", lw=1.0, zorder=3)
        # v3.1.1：tag 只保留 tier，去 ligand 数（信息已下沉到 x 轴副标签）
        tag = f"({tc})"
        dx, dy, ha = label_cfg[g]
        ax.text(x + dx, score + dy, f"$\\it{{{g}}}$ {tag}",
                color="black", fontsize=FS_TICK - 0.5, ha=ha, va="center")

    # 不可成药区（Tdark + Tbio：均无批准药/配体；仅 Tchem 有研究配体）
    ax.axvspan(0.5, 2.5, color="#FFD6D6", alpha=0.4, zorder=1)
    ax.text(1.5, 7.0, "Not directly\ndruggable",
            color="black", fontsize=FS_TICK - 0.5,
            ha="center", va="center", style="italic", weight="bold")

    ax.set_xlim(0.5, 3.7)
    ax.set_ylim(0, 8)
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(["Tdark\n(0 ligand)", "Tbio\n(0 approved)", "Tchem\n(3 ligands)"],
                       fontsize=FS_TICK - 0.5)
    ax.tick_params(axis="x", pad=2)
    ax.set_xlabel("Pharos druggability tier", fontsize=FS_AXIS_LABEL)
    ax.set_ylabel("Intervenability score", fontsize=FS_AXIS_LABEL)
    ax.tick_params(axis="y", labelsize=FS_TICK)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    ax.set_title("c  Direct druggability vs intervenability",
                 fontsize=FS_PANEL_LABEL, loc="left", pad=4, weight="bold")


def panel_d(ax):
    """d. Pharos→表观层出口示意（5-azacytidine 用词=正文 R7）"""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 14)
    ax.axis("off")
    ax.set_xticks([])
    ax.set_yticks([])

    # v3.0 重排：11pt 单行标题超框（"C1GALT1C1 promoter window (reversible)"
    # 542px > box 378px 右溢 168px）且 title/body 垂直仅 8px。改**分行堆叠
    # 引擎**：每层给 title_lines/body_lines 列表，行 pitch 0.80/0.64 unit
    # （行间 ≥6px 视觉分离），内容块在 box 内垂直居中，行高按行数动态算
    # box 高。title 8.5pt / body 7.2pt。三层语义（直接酶抑制红✘ → 转录/
    # 表观层绿✓ → C1GALT1C1 启动子窗深绿✓✓），语义词 5-azacytidine /
    # reversible / [23] 全保留（[23] = x3 docx 线的 Sun 2015，见 v3.6）。
    layers = [
        dict(
            title=["Direct enzyme inhibition"],
            body=["Pharos Tbio/Tdark/Tchem",
                  "0 approved drugs for axis enzymes"],
            bg="#FFE0E0", edge="#D55E00", mark="\u2718",
        ),
        dict(
            title=["Transcriptional /",
                   "epigenetic control"],
            body=["DNMTi (5-aza class)",
                  "EGCG \u00b7 class I HDACi"],
            bg="#D6F0E0", edge="#009E73", mark="\u2713",
        ),
        dict(
            title=["C1GALT1C1 promoter window",
                   "(reversible)"],
            body=["Promoter CpG hypermethylation",
                  "reversible by decitabine [23]"],
            bg="#A8E6CF", edge="#00755E", mark="\u2713\u2713",
        ),
    ]

    PAD_TOP = 0.34
    PAD_BOT = 0.30
    TITLE_LH = 0.86   # title 行 pitch（8.5pt ≈ 18px 字 → 行间 ≥5.5px）
    BODY_LH = 0.72    # body 行 pitch（7.2pt ≈ 15px 字 → 行间 ≥4.8px）
    TITLE_BODY_SEP = 0.16

    def content_pitch(lay):
        return (len(lay["title"]) * TITLE_LH + TITLE_BODY_SEP
                + len(lay["body"]) * BODY_LH)

    def box_height(lay):
        return PAD_TOP + content_pitch(lay) + PAD_BOT

    hs = [box_height(l) for l in layers]        # ≈ 2.86 / 3.66 / 3.66
    GAP = 0.90
    # 从顶向下排布（layers[0]=Direct enzyme inhibition 在最上）：
    # 内容总高 11.98 + 顶边 0.6 + 底边 0.4 → ylim 0-13。
    Y_TOP_CONTENT = 13.4
    bottoms = []
    t = Y_TOP_CONTENT
    for h in hs:
        bottoms.append(t - h)
        t -= h + GAP
    tops = [bottoms[i] + hs[i] for i in range(3)]
    # 校验：bottoms 应 [9.54, 4.98, 0.42] 递减；bottom3=0.42>0

    for i, lay in enumerate(layers):
        yb = bottoms[i]
        yt = tops[i]
        ax.add_patch(plt.Rectangle((0.5, yb), 9.0, yt - yb,
                                   facecolor=lay["bg"], edgecolor=lay["edge"],
                                   lw=1.6, zorder=1))
        # 内容块垂直居中：inner 上界=yt-PAD_TOP、下界=yb+PAD_BOT
        inner_top = yt - PAD_TOP
        inner_bot = yb + PAD_BOT
        nudge = max(0.0, (inner_top - inner_bot - content_pitch(lay)) / 2)
        y_cursor = inner_top - nudge
        for k, tl in enumerate(lay["title"]):
            y_c = y_cursor - (k + 0.5) * TITLE_LH
            ax.text(0.8, y_c, tl, color=lay["edge"],
                    fontsize=8.5, ha="left", va="center",
                    weight="bold", zorder=2)
        y_body0 = y_cursor - len(lay["title"]) * TITLE_LH - TITLE_BODY_SEP
        for k, bl in enumerate(lay["body"]):
            y_c = y_body0 - (k + 0.5) * BODY_LH
            ax.text(0.8, y_c, bl, color="black",
                    fontsize=7.2, ha="left", va="center", zorder=2)
        # 判定标记（右侧；✓✓ 宽 70px≈1.7unit，中心 x=8.55 防右溢框边 9.5）
        ax.text(8.45, (yb + yt) / 2, lay["mark"], color=lay["edge"],
                fontsize=15, ha="center", va="center",
                weight="bold", zorder=2)

    # 箭头（相邻 box 间隙内：上盒底 → 下盒顶）
    for k in range(2):
        y_from = bottoms[k]      # 上盒底
        y_to = tops[k + 1]       # 下盒顶
        ax.annotate("", xy=(5.0, y_to + 0.30), xytext=(5.0, y_from - 0.30),
                    arrowprops=dict(arrowstyle="->", color="black", lw=1.5))

    ax.set_title("d  From Pharos tier to epigenetic window",
                 fontsize=FS_PANEL_LABEL, loc="left", pad=4, weight="bold")


def main():
    # 2x2 网格（加大宽度容纳图例外置；加大高度容纳 panel d 三 box）
    fig = plt.figure(figsize=(mm_to_in(220), mm_to_in(180)), dpi=DPI_DRAFT)
    from matplotlib.gridspec import GridSpec
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1, 1.20],
                  width_ratios=[1, 1], hspace=0.50, wspace=0.40)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    panel_a(ax_a)
    panel_b(ax_b)
    panel_c(ax_c)
    panel_d(ax_d)

    # v3.2 panel a 图例：4 色块 fig-level legend 横向居中于 a 图正下方。
    # a 图位于 gs[0,0]，fig 归一化坐标横向 [0.05, 0.48]，纵向 [0.50, 0.96]
    # （fig 180mm 高，2 行 hspace=0.50）。bbox_to_anchor=(0.265, 0.495)
    # loc='upper center' 会让色块向上展开压 a 图——改 loc='lower center'
    # 让色块向下展开填满 a/c 间隙，归属明确。
    fig.legend(handles=panel_a.legend_handles,
               labels=panel_a.legend_labels,
               loc="lower center",
               bbox_to_anchor=(0.265, 0.485),
               ncol=4,
               frameon=False,
               fontsize=FS_LEGEND - 0.5,
               handlelength=1.6,
               handleheight=1.0,
               columnspacing=1.8)

    out_stem = os.path.join(OUT_DIR, os.environ.get("OUT_STEM", "Fig6_v3.7"))
    save_final(fig, out_stem)
    plt.close(fig)
    print("WROTE", out_stem + ".png")


if __name__ == "__main__":
    main()