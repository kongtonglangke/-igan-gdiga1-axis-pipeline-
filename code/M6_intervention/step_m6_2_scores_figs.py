# -*- coding: utf-8 -*-
"""
M6-2 可干预性评分图
===================
读取 M6_靶点药物证据表.tsv / M6_可干预性评分.tsv，生成三面板评分图：
  A  堆叠水平条：5 基因 × 评分构成（遗传锚点G + 药物手柄D + 功能证据F − 风险扣分P）与总分标注
  B  热图：基因 × 干预窗口/药物机制 的『证据等级』矩阵（0-3 级，与 M6-1 ev_tier 对应）
  C  散点：直接可成药性(Pharos IDG 等级) vs 可干预性总分 —— 展示『轴基因直接不可成药，
     干预杠杆集中于转录/表观窗口（C1GALT1/C1GALT1C1）』
输出：figures/M6_可干预性评分图.png / .pdf
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

OUT = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(OUT, "figures")
os.makedirs(FIGDIR, exist_ok=True)

# ---- 中文字体 ----
for fp in [r"C:/Windows/Fonts/msyh.ttc", r"C:/Windows/Fonts/msyhbd.ttc",
           r"C:/Windows/Fonts/simhei.ttf", r"C:/Windows/Fonts/simsun.ttc"]:
    if os.path.exists(fp):
        font_manager.fontManager.addfont(fp)
plt.rcParams["font.family"] = "Microsoft YaHei"
plt.rcParams["axes.unicode_minus"] = False

# ---- 读取 ----
df_ev = pd.read_csv(os.path.join(OUT, "M6_靶点药物证据表.tsv"), sep="\t")
df_sc = pd.read_csv(os.path.join(OUT, "M6_可干预性评分.tsv"), sep="\t")

gene_order = df_sc["基因"].tolist()  # 总分降序
cmap_comp = {"评分_遗传锚点G(0-3)": "#4C72B0",   # 蓝
             "评分_药物手柄D(0-3)": "#DD8452",   # 橙
             "评分_功能证据F(0-2)": "#55A868"}   # 绿

dfm = df_ev.set_index("基因").loc[gene_order].reset_index()
G = dfm["评分_遗传锚点G(0-3)"].values
D = dfm["评分_药物手柄D(0-3)"].values
F = dfm["评分_功能证据F(0-2)"].values
P = dfm["评分_风险扣分P(0-2)"].values
TOT = G + D + F - P
y = np.arange(len(gene_order))[::-1]

# ---- 面板 B：窗口证据等级矩阵 ----
colsB = ["DNMTi\n(核苷: 5-azaC/DAC)", "DNMTi\n(非核苷: EGCG/RG108)",
         "HDACi\n(class I)", "HDACi\n(泛)", "直接酶抑制剂\n(研究级)"]
mat = {
    "C1GALT1C1": [3, 1, 0, 0, 0],
    "C1GALT1":   [1, 1, 1, 1, 0],
    "GALNT12":   [0, 0, 0, 0, 0],
    "GALNT2":    [0, 0, 0, 0, 2],
    "ST6GALNAC2":[0, 0, 0, 0, 1],
}
matB = np.array([mat[g] for g in gene_order], dtype=float)
# 注：GALNT2 的『直接酶抑制剂』为 GalNAc-T2 研究级选择性抑制剂（JACS Au 2024，IC50≈21 µM），非表观药物
tier_lab = {0: "0 无/不适用", 1: "1 假设/相关性", 2: "2 跨疾病功能或\n研究级工具", 3: "3 IgAN B 细胞\n直接功能"}

pharos_map = {"GALNT2": 2, "C1GALT1": 1, "C1GALT1C1": 1, "ST6GALNAC2": 1, "GALNT12": 0}
pharos_lab = {"GALNT2": "Tchem", "C1GALT1": "Tbio", "C1GALT1C1": "Tbio",
              "ST6GALNAC2": "Tbio", "GALNT12": "Tdark"}

# ================= 绘图 =================
fig = plt.figure(figsize=(17.5, 5.6), constrained_layout=True)
gs = fig.add_gridspec(1, 3, width_ratios=[1.28, 1.18, 1.0], wspace=0.16)

# ---------- A 堆叠条 ----------
axA = fig.add_subplot(gs[0])
left = np.zeros(len(gene_order))
for col, c in cmap_comp.items():
    v = dfm[col].values
    axA.barh(y, v, left=left, color=c, edgecolor="white", linewidth=0.7, label=col.split("(")[0].replace("评分_", ""))
    left += v
axA.barh(y, P, left=left, color="#C44E52", edgecolor="white", linewidth=0.7, hatch="//",
         label="风险扣分 −P")
for i, (tot, pv) in enumerate(zip(TOT, P)):
    yy = y[i]
    axA.text(tot + pv + 0.12, yy, f"{tot:.1f}", va="center", ha="left",
             fontsize=12, fontweight="bold", color="#1a1a1a")
axA.set_yticks(y)
axA.set_yticklabels(gene_order, fontsize=11)
axA.set_xlim(0, 8.5)
axA.set_xlabel("可干预性评分（0–8）", fontsize=10.5)
axA.set_title("A  轴基因可干预性评分（构成分解）", fontsize=11.5, pad=8)
axA.legend(loc="lower right", fontsize=8.6, frameon=True, framealpha=0.95)
axA.grid(axis="x", linestyle=":", alpha=0.35)
axA.set_axisbelow(True)

# ---------- B 热图 ----------
axB = fig.add_subplot(gs[1])
im = axB.imshow(matB, aspect="auto", cmap="YlGnBu", vmin=0, vmax=3)
axB.set_xticks(np.arange(len(colsB)))
axB.set_xticklabels(colsB, fontsize=8.8)
axB.set_yticks(np.arange(len(gene_order)))
axB.set_yticklabels(gene_order, fontsize=11)
for i in range(matB.shape[0]):
    for j in range(matB.shape[1]):
        v = matB[i, j]
        axB.text(j, i, str(int(v)), ha="center", va="center", fontsize=10.5,
                 color="#111" if v in (0, 1) else "white", fontweight="bold")
cb = fig.colorbar(im, ax=axB, fraction=0.046, pad=0.03)
cb.set_ticks([0, 1, 2, 3])
cb.set_ticklabels([tier_lab[0], tier_lab[1], tier_lab[2], tier_lab[3]])
cb.ax.tick_params(labelsize=7.4)
cb.set_label("证据等级", fontsize=9)
axB.set_title("B  基因 × 干预机制 证据等级矩阵", fontsize=11.5, pad=8)
axB.tick_params(top=False)

# ---------- C 散点 ----------
axC = fig.add_subplot(gs[2])
xs = np.array([pharos_map[g] for g in gene_order], dtype=float)
ys = TOT
axC.axvspan(-0.4, 1.4, color="#C44E52", alpha=0.08)
axC.axvline(1.4, color="#C44E52", linestyle="--", linewidth=0.9, alpha=0.7)
colors = [cmap_comp["评分_遗传锚点G(0-3)"] if "G" else "#999" for _ in gene_order]
sc = axC.scatter(xs, ys, s=130, c=colors, edgecolor="#333", linewidth=0.8, zorder=3)
for i, g in enumerate(gene_order):
    dx = 0.0
    if g == "GALNT2":
        dx = 0.0
    axC.annotate(f"{g} ({pharos_lab[g]})", (xs[i], ys[i]),
                 textcoords="offset points", xytext=(10, 6 if i % 2 == 0 else -8),
                 fontsize=9.5, zorder=4)
axC.set_xticks([0, 1, 2])
axC.set_xticklabels(["Tdark\n(0)", "Tbio\n(1)", "Tchem\n(2)"], fontsize=9)
axC.set_xlim(-0.4, 2.8)
axC.set_ylim(1.5, 7.0)
axC.set_xlabel("直接酶位点可成药性（Pharos IDG）", fontsize=10.5)
axC.set_ylabel("可干预性总分（0–8）", fontsize=10.5)
axC.set_title("C  直接可成药性 vs 干预总分", fontsize=11.5, pad=8)
axC.text(0.5, 6.45, "直接不可成药区\n(5/5 基因)", ha="center", va="center",
         fontsize=9, color="#C44E52")
axC.text(2.3, 6.45, "GALNT2 仅\n研究级配体", ha="center", va="center", fontsize=8, color="#555")
axC.grid(linestyle=":", alpha=0.35)

fig.suptitle("", y=1.0)
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(FIGDIR, f"M6_可干预性评分图.{ext}"), dpi=300)
plt.close(fig)
print("figure saved ->", FIGDIR)
