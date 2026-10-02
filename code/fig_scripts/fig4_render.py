"""fig4_render.py — Fig.4 分步中介证据（2026-09-06 步骤 4.3 制作第四图）
作者：AI 王严课题组代办｜依赖：matplotlib/pandas/numpy（fig43 venv）

读：
  阶段3/M4_中介/sigCpG_pos.tsv                      (panel a CpG 位置)
  阶段3/M4_中介/M4_MR_CpG_x_IgAN.tsv                (panel b MR 森林)
  阶段3/M4_中介/M4_coloc_mQTL_x_IgAN.tsv            (panel c coloc IgAN)
  阶段3/M4_中介/M4_coloc_mQTL_x_eQTL.tsv            (panel c coloc eQTL)
  阶段3/M4_中介/M4_方向一致性_lead.tsv              (panel a 方向一致性 beta)

写：
  阶段4/主图/out/Fig4_v3.2.{png,pdf,tiff}（投稿版；4 垂直面板）

版本链：
  v1.0 2026-09-06 步骤 4.3 冻结
  v1.1 2026-09-06 整合审计升版（CpG 去 cg 前缀；y 标签 "·" 分隔；d 精简；
       panel b "11 estimates, all null"、panel a GoDMC mQTL 标签）
  v2.0 2026-09-06 重叠修复整体版（v1.1 件留档）
  v3.0 2026-09-07 用户过审轮（_fig4_check.py 像素级自检驱动，35 处归零）：
      (a) panel a：gene body 带起点 0→1.5 kb（避开启动子盲区框 patch 重叠
          10×92px）；3' 文字 (gene_end+2.2, 0)→(gene_end+0.45, 0.70) 避开
          cg16101574 bar；CpG ID 标签按 rel_kb 排序后**上下两行交替**
          （down -0.62 / up 0.74）消除 cg19603390↔cg19473623 (2.65 kb 间距)
          与 cg17994788↔cg04827551 (7.2 kb) 标签字面重叠；xlim 右扩
          gene_body/1000+10 容最右 ID 完整；
      (b) panel b：**过滤 4 个 IVW P=nan 行** → rows 15→11 与标题
          "11 estimates" 吻合（11 = 5 lead + 5 top + 1 有效 IVW 0.064）；
          P 标签 x 0.55→0.72 + xlim 右扩 1.10（误差棒右端 max 0.501，
          文本 0.72~1.07 完整不裁）；
      (c) panel c：自检已全绿，未改动；
      (d) panel d：box_h 1.30→1.75、gap 0.55→0.60、ylim 0-11→0-13、
          panel_h d 110→150 mm；body 全改单行精简（原 2 行 41px 超 box
          高 39px）；title +0.30→+0.35 / body -0.30→-0.40 分离
          （原互叠 13px）；底部 annotation 文字精简 + y 0.30→0.45。
      自检 _fig4_check.py：text/patch/containment/ticklabels/legend 五件套
      全绿（架构沿 _fig3_check.py，TICK_TOL=3）。
  v3.1 2026-09-07 用户反馈"panel a 仍有几处明显重叠"修复：
      (1) 删 "TSS" / "3'" 文字标记（盲区框 + 3' 端竖线已表达边界）；
      (2) 启动子盲区文字改双行分置盲区框**上下外侧**——旧单行
          (0, 0.68) 同时贴 0.0 刻度 + 与 bar 上行 CpG ID (y=0.74)
          垂直仅 0.06 unit ≈ 2 px 字面紧贴。新版 "TSS ±1.5 kb" y=0.62
          va=bottom 在盲区框上沿 0.50 上方；"(no significant CpG)"
          y=-0.62 va=top 在盲区框下沿 -0.50 下方；
      (3) CpG ID 上行 y=0.74 → 0.92（与盲区文字间距扩到 0.30 unit
          ≈ 10 px），下行同步 -0.62 → -0.78；
      (4) ylim -0.85~0.95 → -1.05~1.20 给双侧 ID 标签 + 盲区上下
          文字留足垂直空间；
      (5) 盲区框上下沿 0.45 → 0.50（容下文字下/上沿）；
      (6) 删除原 (0, -0.22) 副标题 text（与底部 "0" tick 紧邻
          垂直仅 5 px），chr7 坐标 + allele 说明合并进 ylabel 双行；
      (7) 标题缩短 "a Gene-body CpG landscape (none in TSS ±1.5 kb)"
          pad 4→6；
      自检 _fig4_check.py 五件套 TOTAL=0 全绿。
  v3.2 2026-09-07 用户反馈"panel a 还有好几处重叠" + "d 图分色是否
      故意"。_diag4a.py 超严诊断（零容差相交 + NEAR<6px 邻近度）定位
      自检漏掉的视觉粘连根因：
      (1) **AutoLocator tick 越界扩展**：xlim -4~101.7 却生成 -20/120
          tick、ylim -1.05~1.2 生成 -1.5/1.5 tick（nice 边界扩出
          viewLim，越界 tick 画在轴外）。→ set_xticks(0..100 步20)
          + set_yticks(-1.0..1.0 步0.5) 显式固定，tick 数 15→11；
      (2) **盲区文字横向溢出**：x=0 居中的文字半宽 ~8 unit 左溢过
          ax 左缘直逼 ytick 标签列，与 YTICK '-0.5'/'0.5' 横向交叠
          30px + 纵向 gap 0.2/4.4px 视觉粘连 → "TSS ±1.5 kb" 改放
          盲区框右侧外 (x=1.6, ha=left)；"(no significant CpG)" 删
          （标题已含 "none in TSS ±1.5 kb"，属冗余）；
      (3) panel d 颜色经用户确认**保留语义分色**（绿=支持证据 /
          灰=数据缺口 P-only / 朱红=阴性结论），未改。
      诊断后 text 8→6、tick 15→11；_fig4_check.py TOTAL=0 全绿。
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
    FIG_WIDTH_MM, DPI_DRAFT, DPI_TIFF,
    save_final, mm_to_in,
)

PROJECT_ROOT = _paths.LEGACY
F4A_SIGCPG  = PROJECT_ROOT("阶段3/M4_中介/sigCpG_pos.tsv")
F4B_MR      = PROJECT_ROOT("阶段3/M4_中介/M4_MR_CpG_x_IgAN.tsv")
F4C_COLOC1  = PROJECT_ROOT("阶段3/M4_中介/M4_coloc_mQTL_x_IgAN.tsv")
F4C_COLOC2  = PROJECT_ROOT("阶段3/M4_中介/M4_coloc_mQTL_x_eQTL.tsv")
F4A_DIR     = PROJECT_ROOT("阶段3/M4_中介/M4_方向一致性_lead.tsv")
OUT_DIR     = PROJECT_ROOT("阶段4/主图/out")
os.makedirs(OUT_DIR, exist_ok=True)


def _fmt_p(p):
    try:
        p = float(p)
    except Exception:
        return "NA"
    if p < 1e-4:
        return f"{p:.1e}".replace("e-0", "e-")
    if p < 0.001:
        return f"{p:.1e}"
    if p < 0.01:
        return f"{p:.4f}"
    return f"{p:.3f}"


def new_fig_four_panels(panel_h_mm=(58, 58, 58, 80), w_mm=FIG_WIDTH_MM):
    """4 垂直面板 figure。"""
    total_h_mm = sum(panel_h_mm) + 14
    fig, axes = plt.subplots(
        4, 1,
        figsize=(mm_to_in(w_mm), mm_to_in(total_h_mm)),
        dpi=DPI_DRAFT,
        constrained_layout=True,
        gridspec_kw={"height_ratios": panel_h_mm, "hspace": 0.45},
    )
    return fig, axes


def panel_a(ax):
    """a. C1GALT1 gene body CpG landscape（v3.1 2026-09-07 19:50）

    v3.6 2026-09-11 图文一致性修复（阶段4 ↔ 阶段4.5 双向审查，初稿→包方向）：
      • panel d 步骤④ accession "GCST90018888" → **"GCST90011884"**——Han-Chinese
        定量 Gd-IgA1 GWAS 的正确编号（Methods / Table 1 一致；GCST90018866 是
        Sakaue 2021 IgAN meta），原号为全文任何位置均不存在的错号；
      • panel c 标题 "all <0.014" → **"all ≤0.014"**——与图注 "PPH4 ≤ 0.014"
        及 panel d 步骤⑤ "≤0.014" 统一（实测最大值 0.01383，两式均成立）；
      • panel d 步骤① "P≤1e-321" → **"P ≈ 1e-321"**——底层数据为 1.34e-321，
        "≤1e-321" 方向不成立。
    未改动任何数据、版式几何或结论。"""
    df = pd.read_csv(F4A_SIGCPG, sep="\t")
    TSS = 7196565
    gene_end = 7288282
    gene_body = gene_end - TSS

    df["rel_kb"] = (df["pos37"] - TSS) / 1000.0
    direction = pd.read_csv(F4A_DIR, sep="\t")
    dir_r13 = direction[(direction["lead"] == "rs13226913") &
                       (direction["context"] == "GTEx_v10_blood")]
    beta_dict = dict(zip(dir_r13["cpg"], dir_r13["mqtl_beta_on_eQTLalt"]))

    # 基因体背景：从 TSS+1.5 kb 起（0~1.5 为启动子盲区，由下方 hatch 框
    # 标注，避免基因带与盲区框 patch 重叠 10×92 px）。
    ax.add_patch(plt.Rectangle((1.5, -0.45), gene_body/1000.0 - 1.5, 0.9,
                               facecolor="#E5F1F7", edgecolor="#A8C8DC",
                               lw=0.8, zorder=1))
    # TSS 竖线（无 "TSS" 文字，避免与 0.0 刻度/盲区文字/标题紧贴）
    ax.axvline(0, color="black", lw=1.5, zorder=3)
    # 3' 竖线（无 "3'" 文字；cg16101574 bar 顶 0.43 上方已有
    # CpG ID 上行标签位于 y=0.92，"3'" 文字属冗余）
    ax.axvline(gene_body/1000.0, color="black", lw=1.5, zorder=3)
    # 启动子盲区（TSS ±1.5 kb）
    # v3.1：盲区文字双行分置框上下外侧 (x=0 居中, y=±0.62)
    # v3.2：文字全改放盲区框**右侧外** (x=1.6 ha=left 单行)——x=0
    #   居中时文字半宽 ~8 unit 左溢过 ax 左缘直逼 ytick 标签列（与
    #   YTICK '-0.5'/'0.5' 横向交叠 30px + 纵向 gap <5px，_diag4a.py
    #   实测 NEAR 0.2px/4.4px 视觉粘连）。"(no significant CpG)" 删
    #   （标题已含 "none in TSS ±1.5 kb"），只留 "TSS ±1.5 kb" 归属
    #   标注从 x=1.6 向右展开，不碰任何 ytick/ID 标签。
    ax.add_patch(plt.Rectangle((-1.5, -0.50), 3.0, 1.0,
                               facecolor="none", edgecolor="black",
                               hatch="\\\\", lw=0.8, zorder=2))
    ax.text(1.6, 0.55, "TSS ±1.5 kb", color="black",
            fontsize=FS_TICK - 0.8, ha="left", va="center",
            style="italic")

    # 5 CpG 标记 + ID 标签：v3.0 按 rel_kb 排序后 ID 标签**上下交替两行**
    # （down y=-0.62 / up y=0.74）。根因：ID 8 字符标签像素宽 ~80 px >
    # 相邻 CpG 间距（cg19603390↔cg19473623 仅 2.65 kb ≈ 17 px →
    # 重叠 62×16px；cg17994788↔cg04827551 7.2 kb → 重叠 31×16px）。
    # 交替后同层最近对间距 cg19603390→cg17994788 39 kb / up 层 43.9 kb
    # 均 > 标签宽对应 kb，两两无叠。up 层 y=0.74 高于 bar 顶
    # (β max 0.43) 与 TSS±1.5 文字 y=0.68 区无叠（x 距 >25 kb）。
    # v3.3：up 层 y=0.92→0.78（避开 cg19473623 在 β=0.43 时柱顶 y≈0.68 + 字高
    # 约 0.10 → 字底 y≈0.58；y=0.78 时字底 ≈0.68 距 β=0.43 顶 0.10 unit 余量）；
    # down 层同步 -0.78→-0.85；ylim 扩 -1.10~1.18 给下层 ID 留 0.10 unit padding
    df = df.sort_values("rel_kb").reset_index(drop=True)
    for i, row in df.iterrows():
        x = row["rel_kb"]
        cpg = row["cpg"]
        short = cpg[2:]  # 去 "cg"
        beta = beta_dict.get(cpg, 0.0)
        color = OKABE_ITO["vermillion"] if beta >= 0 else OKABE_ITO["blue"]
        ax.bar(x, beta, width=0.5, bottom=0, color=color,
               edgecolor="black", lw=0.4, zorder=4)
        # v3.3：up 层 0.92→0.78 + down 层 -0.78→-0.85（更贴近数据区，避开 ax 边界）
        if i % 2 == 0:
            ax.text(x, -0.85, short, color="black", fontsize=FS_TICK - 0.8,
                    ha="center", va="top", rotation=0)
        else:
            ax.text(x, 0.78, short, color="black", fontsize=FS_TICK - 0.8,
                    ha="center", va="bottom", rotation=0)

    ax.set_xlim(-4, gene_body/1000.0 + 10)
    ax.set_ylim(-1.10, 1.18)
    # v3.2：显式固定刻度——AutoLocator 会把 nice 边界扩展到 viewLim
    #   之外（xlim -4~101.7 却生成 -20/120 tick，ylim -1.05~1.2 生成
    #   -1.5/1.5 tick，越界刻度画在轴外侧视觉杂乱 + 与轴角文字近贴）。
    ax.set_xticks(np.arange(0, 101, 20))     # 0,20,40,60,80,100
    ax.set_xticklabels(["0", "20", "40", "60", "80", "100"],
                       fontsize=FS_TICK)
    ax.set_yticks([-1.0, -0.5, 0.0, 0.5, 1.0])
    ax.set_yticklabels(["-1.0", "-0.5", "0.0", "0.5", "1.0"],
                       fontsize=FS_TICK)
    ax.set_ylabel("cis-mQTL β (chr7:7.20–7.29 Mb; GoDMC whole blood;\n"
                  "sign per eQTL-raising allele of rs13226913)",
                  fontsize=FS_AXIS_LABEL, labelpad=8)
    ax.tick_params(axis="both", labelsize=FS_TICK, pad=2)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    ax.axhline(0, color="black", lw=0.5, zorder=2)
    ax.set_title("a  Gene-body CpG landscape (none in TSS ±1.5 kb)",
                 fontsize=FS_PANEL_LABEL, loc="left", pad=6, weight="bold")
    # v3.1 删除原 (0, -0.22) 副标题 ax.text ——与底部 "0" tick 紧邻
    # （垂直 5 px），且新版 ylabel 已合并 chr7 坐标 + allele 说明，
    # 副标题属冗余。


def panel_b(ax):
    """b. CpG → IgAN MR forest (all null)"""
    df = pd.read_csv(F4B_MR, sep="\t")
    rows = []
    for cpg in df["cpg"].unique():
        sub = df[df["cpg"] == cpg]
        r_lead = sub[sub["method"] == "Wald_lead"]
        r_topcis = sub[sub["method"] == "Wald_topcis"]
        r_ivw = sub[sub["method"].str.startswith("IVW")]
        # v3.0：仅 append P 有限值行。根因：4/5 CpG 的 IVW P=nan 也被
        # 当作 estimate 绘制 → rows=15，P 标签行距过挤互相 1px 重叠，
        # 且与标题 "11 estimates" 矛盾（11 = 5 lead + 5 top + 1 有效 IVW，
        # cg04827551 P=0.064）。P nan 行跳过，行距恢复，标题数字吻合。
        if not r_lead.empty:
            r1 = r_lead[r_lead["instrument"] == "rs13226913"]
            if not r1.empty:
                p_lead = float(r1.iloc[0]["mr_p"])
                if np.isfinite(p_lead):
                    rows.append((cpg, "lead", float(r1.iloc[0]["mr_beta"]),
                                 float(r1.iloc[0]["mr_se"]), p_lead))
        if not r_topcis.empty:
            p_top = float(r_topcis.iloc[0]["mr_p"])
            if np.isfinite(p_top):
                rows.append((cpg, "top", float(r_topcis.iloc[0]["mr_beta"]),
                             float(r_topcis.iloc[0]["mr_se"]), p_top))
        if not r_ivw.empty:
            p_ivw = float(r_ivw.iloc[0]["mr_p"])
            if np.isfinite(p_ivw):
                rows.append((cpg, "ivw", float(r_ivw.iloc[0]["mr_beta"]),
                             float(r_ivw.iloc[0]["mr_se"]), p_ivw))

    n = len(rows)
    y = np.arange(n)[::-1]
    for i, (cpg, label, beta, se, p) in enumerate(rows):
        ax.errorbar(beta, y[i], xerr=1.96*se, fmt="o", color=OKABE_ITO["blue"],
                    capsize=3, markersize=5, lw=1.2, zorder=3)
        # P 标签：v3.0 统一右列 x=0.72 起（数据误差棒右端 max ≈0.50，
        # cg17994788 top 0.501），xlim 右扩 1.10 完整容纳 "P=0.819" 等
        # 8 字符文本（~0.35 unit 宽 → 0.72~1.07）。旧 x=0.55 与 xlim 0.78
        # 右缘裁切文本 + 与 0.501 cap 间隙过窄。
        ax.text(0.72, y[i], f"P={_fmt_p(p)}", color=OKABE_ITO["blue"],
                fontsize=FS_TICK - 1, va="center", ha="left")

    ax.axvline(0, color="black", lw=0.8, ls="--", zorder=2)
    ax.set_xlim(-0.4, 1.10)
    ax.set_xticks([-0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xticklabels(["-0.4", "-0.2", "0.0", "0.2", "0.4", "0.6", "0.8", "1.0"],
                       fontsize=FS_TICK)
    ax.set_yticks(y)
    # Y 标签：精简为单行（CpG 短 ID + 方法 括号；，避免 2 行重叠）
    yt = [f"$\\it{{{cpg[2:]}}}$ ({lbl})" for cpg, lbl, *_ in rows]
    ax.set_yticklabels(yt, fontsize=FS_TICK - 1.5)
    ax.tick_params(axis="y", pad=2)
    ax.set_xlabel("MR β (CpG methylation → IgAN; Wald or IVW)",
                  fontsize=FS_AXIS_LABEL)
    ax.tick_params(axis="x", labelsize=FS_TICK)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    ax.set_title("b  CpG → IgAN Mendelian randomization (11 estimates, all null)",
                 fontsize=FS_PANEL_LABEL, loc="left", pad=4, weight="bold")


def panel_c(ax):
    """c. Coloc PPH4: mQTL × IgAN/eQTL_GTEx/eQTL_B_naive"""
    igan = pd.read_csv(F4C_COLOC1, sep="\t")
    eqtl = pd.read_csv(F4C_COLOC2, sep="\t")
    cpgs = ["cg19603390", "cg19473623", "cg17994788", "cg04827551", "cg16101574"]

    pph_igan = [float(igan[igan["cpg"] == c].iloc[0]["PPH4"]) for c in cpgs]
    eqtl_blood = eqtl[eqtl["context"] == "GTEx_v10_blood"]
    eqtl_kidney = eqtl[eqtl["context"] == "GTEx_v10_kidney_cortex"]
    eqtl_naive = eqtl[eqtl["context"] == "OneK1K_B_naive"]
    pph_gblood = [float(eqtl_blood[eqtl_blood["cpg"] == c].iloc[0]["PPH4"]) for c in cpgs]
    pph_gkidney = [float(eqtl_kidney[eqtl_kidney["cpg"] == c].iloc[0]["PPH4"]) for c in cpgs]
    pph_naive = [float(eqtl_naive[eqtl_naive["cpg"] == c].iloc[0]["PPH4"]) for c in cpgs]

    x = np.arange(len(cpgs))
    width = 0.20

    ax.bar(x - 1.5*width, pph_igan, width=width, color=OKABE_ITO["vermillion"],
           edgecolor="black", lw=0.4, label="mQTL×IgAN", zorder=3)
    ax.bar(x - 0.5*width, pph_gblood, width=width, color=OKABE_ITO["blue"],
           edgecolor="black", lw=0.4, label="mQTL×GTEx blood", zorder=3)
    ax.bar(x + 0.5*width, pph_gkidney, width=width, color=OKABE_ITO["sky_blue"],
           edgecolor="black", lw=0.4, label="mQTL×GTEx kidney", zorder=3)
    ax.bar(x + 1.5*width, pph_naive, width=width, color=OKABE_ITO["orange"],
           edgecolor="black", lw=0.4, label="mQTL×B naive", zorder=3)

    ax.axhline(0.05, color=C_GW_LINE, lw=1.0, ls="--", zorder=2)
    ax.text(0.01, 0.12, "PPH4 = 0.05", color=C_GW_LINE,
            fontsize=FS_TICK - 0.5, ha="left", va="bottom",
            transform=ax.transAxes)

    ax.set_xticks(x)
    # X 标签去 "cg" 前缀缩小 + 旋转 35° + ha='right'（避免水平叠加）
    ax.set_xticklabels([c[2:] for c in cpgs], fontsize=FS_TICK - 0.8,
                       rotation=35, ha="right")
    ax.tick_params(axis="x", pad=2)
    ax.set_ylabel("PPH4 (coloc posterior)", fontsize=FS_AXIS_LABEL)
    ax.set_ylim(0, 1.0)
    ax.tick_params(axis="y", labelsize=FS_TICK)
    # v3.3：图例改 ax 内右下空白区（避开 cg17994788/cg04827551 高柱
    # [x=2,3]; 此前 ax 右侧外 (1.02,1.00) 受 constrained_layout 限制，
    # 图例被推到 ax 内左上挡柱）。ncol=1 紧凑 + frameon 白底遮柱子。
    ax.legend(loc="upper left", bbox_to_anchor=(0.02, 0.98),
              fontsize=FS_LEGEND - 0.3, frameon=True, framealpha=0.95,
              ncol=1, borderaxespad=0.4)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    ax.set_title("c  Coloc PPH4: mQTL→IgAN all ≤0.014; mQTL→B-naive eQTL ≥0.79",
                 fontsize=FS_PANEL_LABEL, loc="left", pad=4, weight="bold")


def panel_d(ax):
    """d. Stepwise mediation framework with proportion gap"""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 13)
    ax.axis("off")

    # v3.0 重排：box_h 1.30→1.75 / gap 0.55→0.60（原 box 高 ~33 px 装不下
    # title 11pt(23px)+body 2 行(41px) → 文字上下溢出 5~14 px 且 title/body
    # 在 box 内互叠 13px）。body 全部改**单行精简**（原 2 行文字让 box 内
    # 垂直空间必然不够）。ylim 0-11→0-13 + panel_h d 110→150 mm 提供
    # 垂直余量；title y_pos+0.35 / body y_pos-0.40 分离。
    # 精简正文（单行，最长 ~45 字符，8pt 单行 bbox < box 宽 9.4 unit）：
    boxes = [
        ("① Lead SNP → gene-body CpG (mQTL)",
         "GoDMC n=25–28k; P ≈ 1e-321 (5 CpG)",
         OKABE_ITO["bluish_grn"]),
        ("② CpG → C1GALT1 expression",
         "B naive PPH4=0.79–0.82 (2 CpG)",
         OKABE_ITO["bluish_grn"]),
        ("③ Lead SNP → C1GALT1 (cis-eQTL)",
         "β>0; rs13226913 B naive P=9.4e-6",
         OKABE_ITO["bluish_grn"]),
        ("④ C1GALT1 → Gd-IgA1 (quantitative)",
         "GCST90011884; P-only (no β)*",
         OKABE_ITO["gray"]),
        ("⑤ CpG → IgAN (Mendelian randomization)",
         "MR all null; coloc PPH4≤0.014",
         OKABE_ITO["vermillion"]),
    ]

    box_h = 1.75
    gap = 0.60
    y_top = 12.3
    for i, (title, body, color) in enumerate(boxes):
        y_pos = y_top - i * (box_h + gap) - box_h/2
        rect = plt.Rectangle((0.3, y_pos - box_h/2), 9.4, box_h,
                              facecolor="white", edgecolor=color, lw=2.0)
        ax.add_patch(rect)
        # 标题：box 中心偏上（+0.35），避开下方 body
        ax.text(5.5, y_pos + 0.35, title, color=color,
                fontsize=FS_TITLE, ha="center", va="center", weight="bold")
        # 正文：box 中心偏下（-0.40），单行
        ax.text(5.5, y_pos - 0.40, body, color="black",
                fontsize=FS_TICK - 0.5, ha="center", va="center",
                style="italic")
        if i < len(boxes) - 1:
            # 箭头：连接相邻 box 的下沿→上沿（box 底 = y_pos-box_h/2；
            # 下 box 顶 = y_pos-box_h/2-gap，gap-0.05 微缩不触边）
            arrow_y_from = y_pos - box_h/2
            arrow_y_to = y_pos - box_h/2 - gap + 0.05
            ax.annotate("", xy=(5.5, arrow_y_to), xytext=(5.5, arrow_y_from),
                        arrowprops=dict(arrowstyle="->", color="black", lw=1.2))

    # 底部 annotation：放 box 4 底 (1.15) 之下，y=0.45 中心（bbox 外框
    # 顶 <1.15 不触框）；文字 2 行精简
    ax.text(5.0, 0.45,
            "Mediation proportion: not estimable\n(P-only at step ④)",
            color=C_GW_LINE, fontsize=FS_TICK - 0.3,
            ha="center", va="center", weight="bold",
            bbox=dict(boxstyle="round,pad=0.35", facecolor="white",
                      edgecolor=C_GW_LINE, lw=1.0))

    ax.set_title("d  Stepwise mediation framework (proportion gap highlighted)",
                 fontsize=FS_PANEL_LABEL, loc="left", pad=4, weight="bold")


def main():
    # v3.0：panel d 110→150 mm（box 加高后垂直余量）；fig 总高 364 mm
    fig, axes = new_fig_four_panels(panel_h_mm=(60, 70, 70, 150))
    panel_a(axes[0])
    panel_b(axes[1])
    panel_c(axes[2])
    panel_d(axes[3])

    out_stem = os.path.join(OUT_DIR, os.environ.get("OUT_STEM", "Fig4_v3.2"))
    save_final(fig, out_stem)
    plt.close(fig)
    print("WROTE", out_stem + ".png")


if __name__ == "__main__":
    main()