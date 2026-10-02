#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""M4-3 整合中介结果表 + 分步中介证据图"""
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
from matplotlib.patches import FancyBboxPatch

plt.rcParams["font.family"] = "Microsoft YaHei"
plt.rcParams["axes.unicode_minus"] = False

OUT = _paths.at("阶段3/M4_中介")
FIG = OUT("figures")
os.makedirs(FIG, exist_ok=True)

qc = pd.read_csv(OUT("M4_输入层_QC.tsv"), sep="\t")
mr = pd.read_csv(OUT("M4_MR_CpG_x_IgAN.tsv"), sep="\t")
col_i = pd.read_csv(OUT("M4_coloc_mQTL_x_IgAN.tsv"), sep="\t")
col_e = pd.read_csv(OUT("M4_coloc_mQTL_x_eQTL.tsv"), sep="\t")
dirs = pd.read_csv(OUT("M4_方向一致性_lead.tsv"), sep="\t")

# ---------- 1) 中介结果表 ----------
rows = []
for _, r in qc.iterrows():
    cpg = r["cpg"]
    # MR: topcis & lead Wald
    m = mr[mr["cpg"] == cpg]
    def pick(method):
        s = m[m["method"].str.startswith(method)]
        return s.iloc[0] if len(s) else None
    t = pick("Wald_topcis"); l1 = pick("Wald_lead_rs13226913"); l2 = pick("Wald_lead_rs10238682")
    ci = col_i[col_i["cpg"] == cpg].iloc[0]
    # eQTL coloc 主导假设 (per tissue) 摘要
    ce = col_e[col_e["cpg"] == cpg]
    def dom(x):
        pps = {"PPH0": x["PPH0"], "PPH1": x["PPH1"], "PPH2": x["PPH2"],
               "PPH3": x["PPH3"], "PPH4": x["PPH4"]}
        k = max(pps, key=pps.get)
        return "%s=%.3f" % (k, pps[k])
    ce_sum = "; ".join("%s:%s" % (x["context"].replace("GTEx_v10_", "GTEx-").replace("OneK1K_", "1K1K-"),
                                  dom(x)) for _, x in ce.iterrows())
    # 方向一致性: 显著 eQTL 语境下同向比例
    dd = dirs[(dirs["cpg"] == cpg) & (dirs["eqtl_P"] < 0.05)]
    n_same = int(dd["sign_consistent"].sum())
    dsum = "rs13226913:%d/%d; " % (n_same, len(dd)) if len(dd) else ""
    dd2 = dirs[(dirs["cpg"] == cpg) & (dirs["lead"] == "rs10238682") & (dirs["eqtl_P"] < 0.05)]
    dsum += "rs10238682:%d/%d" % (int(dd2["sign_consistent"].sum()), len(dd2)) if len(dd2) else ""
    rows.append(dict(
        cpg=cpg, pos37=r["pos37"], feature=r["feature"], dist_TSS_bp=r["dist_TSS"],
        samplesize=r["sample_n"], n_sig_cis_IV=r["n_sig_cis_IV"], n_clumped_IV=r["n_clumped_IV"],
        top_cis_rsid=r["top_cis_rsid"], top_cis_P=r["top_cis_P"], top_cis_beta=r["top_cis_beta"],
        lead_mQTL="rs13226913+rs10238682",
        # CpG→IgAN MR (Wald top cis)
        MR_CpGtoIgAN_beta=np.nan if t is None else t["mr_beta"],
        MR_CpGtoIgAN_se=np.nan if t is None else t["mr_se"],
        MR_CpGtoIgAN_P=np.nan if t is None else t["mr_p"],
        MR_lead_rs13226913_P=np.nan if l1 is None else l1["mr_p"],
        MR_lead_rs10238682_P=np.nan if l2 is None else l2["mr_p"],
        coloc_mQTLxIgAN_PPH1=ci["PPH1"], coloc_mQTLxIgAN_PPH4=ci["PPH4"],
        IgAN_region_topP=ci["igan_top_P"], n_shared_IgAN=ci["n_shared"],
        coloc_mQTLxeQTL_主导=ce_sum,
        方向一致性_显著eQTL语境=dsum))
res = pd.DataFrame(rows)
res.to_csv(OUT("中介结果表.tsv"), sep="\t", index=False)
print("saved 中介结果表.tsv (%d rows)" % len(res))

# ---------- 2) 图: 分步中介证据 ----------
leads_pos = {"rs13226913": 7_246_846, "rs10238682": 7_255_017}
leads_color = {"rs13226913": "#C0392B", "rs10238682": "#2874A6"}
cpg_sorted = ["cg19603390", "cg19473623", "cg17994788", "cg04827551", "cg16101574"]
raw = pd.read_csv(OUT("lead_x_mQTL_sigCpG_raw.tsv"), sep="\t")

fig = plt.figure(figsize=(13.5, 9.5))
gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], hspace=0.42, wspace=0.30,
                      left=0.09, right=0.97, top=0.93, bottom=0.09)

# --- A: 基因体示意 + lead mQTL 效应 (全血) ---
ax = fig.add_subplot(gs[0, 0])
gs0, ge0, tss0 = 7_196_565, 7_288_282, 7_196_565
ax.add_patch(plt.Rectangle((gs0, -0.28), ge0 - gs0, 0.56, facecolor="#EAF2F8",
                           edgecolor="#2C3E50", lw=1.2, zorder=1))
ax.annotate("", xy=(tss0, 0.42), xytext=(tss0, 0.06),
            arrowprops=dict(arrowstyle="-|>", color="#2C3E50", lw=1.4))
ax.text(tss0 - 7000, 0.50, "TSS\n(+链)", ha="center", fontsize=8)
xr = (gs0 - 30000, ge0 + 40000)
for rs, pos in leads_pos.items():
    ax.plot([pos, pos], [-0.28, 0.30], ls="--", lw=0.9, color=leads_color[rs], alpha=0.6)
    ax.text(pos, 0.44, rs, rotation=90, fontsize=7.5, color=leads_color[rs], ha="center", va="bottom")
ymax = 0
posmap = {"cg19603390": 7_222_222, "cg19473623": 7_224_869, "cg17994788": 7_261_594,
          "cg04827551": 7_268_805, "cg16101574": 7_291_514}
betamap = {}
for _, r in raw.iterrows():
    betamap.setdefault(r["lead_rsid"], {})[r["cpg"]] = float(r["beta_a1"])
for cpg in cpg_sorted:
    for rs in leads_pos:
        b = betamap.get(rs, {}).get(cpg)
        if b is None:
            continue
        ax.bar(posmap[cpg], b, width=1800, color=leads_color[rs], alpha=0.85,
               edgecolor="white", lw=0.4)
        ymax = max(ymax, abs(b))
ax.axhline(0, color="grey", lw=0.6)
ax.set_ylim(-0.5, 1.05)
ax.set_xlim(*xr)
ax.set_xticks([gs0, tss0, 7_222_222, 7_246_846, 7_255_017, 7_291_514, ge0])
ax.set_xticklabels(["7.197", "TSS", "7.222", "7.247", "7.255", "7.292", "7.288"], fontsize=7)
for lbl in ax.get_xticklabels():
    lbl.set_rotation(0)
ax.set_ylabel("全血 mQTL β (每等位, β-值单位)", fontsize=9)
ax.set_title("A. C1GALT1 位点：lead SNP → 基因体/3'侧 CpG 甲基化（GoDMC 全血, n=25–28k）",
             fontsize=10.5, loc="left")
ax.text((gs0 + ge0) / 2, 0.70, "C1GALT1  基因体 (hg19 chr7:7,196,565–7,288,282)", ha="center", fontsize=9.5,
        color="#1B2631")
ax.legend(handles=[plt.Line2D([0], [0], color=leads_color[k], lw=4, label=k) for k in leads_color],
          fontsize=8.5, loc="lower right", frameon=True, framealpha=0.9)
ax.annotate("启动子窗\nTSS±1.5kb\n(无 mQTL 靶点)", xy=(7_197_400, 0.85), fontsize=8, color="#7F8C8D", ha="center",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85))

# --- B: CpG→IgAN MR 森林图 (Wald top-cis) ---
ax = fig.add_subplot(gs[0, 1])
mr2 = mr[mr["method"].str.startswith("Wald_topcis")].set_index("cpg").reindex(cpg_sorted)
ys = np.arange(len(cpg_sorted))[::-1]
for y, cpg in zip(ys, cpg_sorted):
    r_ = mr2.loc[cpg]
    b, s, p = r_["mr_beta"], r_["mr_se"], r_["mr_p"]
    ax.errorbar(b, y, xerr=1.96 * s, fmt="o", color="#2874A6", ms=6, capsize=4, lw=1.3)
    ax.text(0.42, y, "β=%.3f\nP=%.2f" % (b, p), fontsize=8, va="center")
ax.axvline(0, color="grey", lw=0.8)
ax.set_yticks(ys); ax.set_yticklabels(cpg_sorted, fontsize=9)
ax.set_xlim(-0.55, 0.55)
ax.set_xlabel("CpG 甲基化 → IgAN 因果效应 (Wald, 每 β-值单位, 全血 cis-IV)", fontsize=9)
ax.set_title("B. 甲基化→IgAN 因果效应（GCST90018866）：全阴性", fontsize=10.5, loc="left")
ax.text(0.02, 0.03, "与主证据3 '机制轴不驱动 IgAN 易感' 一致\n(甲基化侧延伸验证)", transform=ax.transAxes,
        fontsize=8.5, color="#C0392B")

# --- C: coloc 汇总 (mQTL × IgAN / × eQTL) ---
ax = fig.add_subplot(gs[1, 0])
xpos = np.arange(len(cpg_sorted))
ci_idx = col_i.set_index("cpg").reindex(cpg_sorted)
w = 0.38
pph1 = ci_idx["PPH1"].fillna(0).values
pph4i = ci_idx["PPH4"].fillna(0).values
b1 = ax.bar(xpos - w / 2, pph1, w, color="#E74C3C", alpha=0.85, label="mQTL×IgAN  PPH1(仅甲基化信号)")
b2 = ax.bar(xpos + w / 2, pph4i, w, color="#F5B041", alpha=0.85, label="mQTL×IgAN  PPH4(共享因果)")
# eQTL coloc: PPH3 最高语境比例 → 次轴简化, 打印文本
ce = col_e.set_index(["cpg", "context"]).reindex(
    pd.MultiIndex.from_product([cpg_sorted, col_e["context"].unique()]))
txt = []
for i, cpg in enumerate(cpg_sorted):
    sub = col_e[col_e["cpg"] == cpg]
    mx = sub.loc[sub[["PPH0", "PPH1", "PPH2", "PPH3", "PPH4"]].max(axis=1).idxmax()]
    domk = max(["PPH0", "PPH1", "PPH2", "PPH3", "PPH4"], key=lambda k: mx[k])
    txt.append("%s→%s" % (mx["context"].replace("GTEx_v10_", "GTEx-").replace("OneK1K_", "1K1K-"), domk))
for i, t in enumerate(txt):
    ax.text(i, 1.045, t, rotation=70, ha="right", va="bottom", fontsize=6.8, color="#1A5276")
ax.set_xticks(xpos); ax.set_xticklabels(cpg_sorted, fontsize=9)
ax.set_ylim(0, 1.18)
ax.set_ylabel("coloc 后验概率", fontsize=9)
ax.set_title("C. 共定位：甲基化与 IgAN 无共享信号；与 C1GALT1 表达多为‘区域双信号(PPH3)’",
             fontsize=10.5, loc="left", pad=10)
ax.legend(fontsize=8, frameon=True, framealpha=0.9, loc="lower right")
ax.text(0, -0.20, "注: GoDMC assoc_meta 仅收录显著 SNP–CpG 对, coloc 结果作单信号指示解读(见方法学说明)",
        transform=ax.transAxes, fontsize=7.5, color="#7F8C8D")

# --- D: 分步中介证据链 + Gd-IgA1 缺口 ---
ax = fig.add_subplot(gs[1, 1])
ax.set_xlim(0, 1); ax.set_ylim(0, 1.3)
ax.axis("off")
steps = [
    ("① lead SNP → 基因体 CpG 甲基化", "证据充分 (P≤1e-321, n=25k)", "#27AE60", True),
    ("② 甲基化 → C1GALT1 表达耦合", "部分证据 (方向同向为主, coloc 双信号)", "#F39C12", True),
    ("③ lead SNP → C1GALT1 表达", "证据充分 (血/肾/B细胞, P≤1e-11)", "#27AE60", True),
    ("④ 表达 → Gd-IgA1 (中介比例)", "数据缺口: GCST90011884 仅 P 值", "#E74C3C", False),
    ("⑤ 甲基化 → IgAN 疾病", "阴性 (MR+coloc, 主证据3 延伸)", "#E74C3C", True),
]
for i, (t, s, c, ok) in enumerate(steps):
    y = 1.05 - i * 0.22
    ax.add_patch(FancyBboxPatch((0.04, y - 0.085), 0.92, 0.17,
                 boxstyle="round,pad=0.01", fc="white", ec=c, lw=1.6))
    ax.text(0.07, y + 0.02, t, fontsize=9.3, va="center", color="#1B2631")
    ax.text(0.07, y - 0.05, s, fontsize=7.7, va="center", color=c)
ax.text(0.5, 1.22, "D. 分步中介证据链与正式中介比例的缺口", ha="center", fontsize=11, fontweight="bold")
ax.text(0.5, -0.05, "完成④ 的路径: Kiryluk 2017 Gd-IgA1 GWAS (dbGaP phs000431, 需授权) 或 Wang 2021 作者索取全量效应量",
        ha="center", fontsize=8.2, color="#7F8C8D")
fig.suptitle("M4 甲基化中介分析（公开数据可计算范围）— 阶段3 步骤3.1", fontsize=13, y=0.985)
fig.savefig(os.path.join(FIG, "M4_分步中介证据图.png"), dpi=200)
fig.savefig(os.path.join(FIG, "M4_分步中介证据图.pdf"))
print("saved M4_分步中介证据图.png/pdf")
