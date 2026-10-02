# -*- coding: utf-8 -*-
"""M5 东亚异质性 · gnomAD v4 EAS/EUR 频率回填与对比图
数据源：gnomAD GraphQL API (dataset=gnomad_r4, GRCh38), genome 数据集;
EUR 以 NFE (Non-Finnish European) 为代表; X 位点按 gnomAD X 非 PAR 口径(含男半合子)。
跨祖先 meta EAF 来自 GCST90018866 (IgAN 跨祖先 meta)。
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False

ROWS = [
    # rsid gene gnomad_str alt EAS_AF EAS_AC EAS_AN EUR_AF EUR_AC EUR_AN overall_AF igan_meta_EAF wang_P igan_meta_P
    ["rs13226913", "C1GALT1",   "7-7207215-T-C",  "C", 0.92623, 4796, 5178, 0.42047, 28556, 67914, 0.53297, 0.57614, 0.01088, 0.2842],
    ["rs10238682", "C1GALT1",   "7-7215386-A-G",  "G", 0.52226, 2698, 5166, 0.20354, 13839, 67990, 0.30101, 0.28032, 1.203e-09, 0.222],
    ["rs7856182",  "GALNT12",   "9-98871030-C-T", "T", 0.02530, 131,  5178, 0.15404, 10472, 67982, 0.19338, 0.11379, 2.379e-09, 0.1004],
    ["rs5910940",  "C1GALT1C1", "X-120624955-G-A","A", 0.54648, 1934, 3539, 0.51109, 27077, 52979, 0.55382, float("nan"), 0.004042, float("nan")],
]

df = pd.DataFrame(ROWS, columns=["lead_rsid","lead_gene","gnomad_query","alt",
                                 "gnomAD_EAS_AF","EAS_AC","EAS_AN",
                                 "gnomAD_EUR_NFE_AF","EUR_AC","EUR_AN",
                                 "gnomAD_overall_AF","IgAN_meta_EAF","Wang2021_Han_P","IgAN_meta_P"])
df.to_csv("M5_gnomAD_EASEUR_频率_20260903.tsv", sep="\t", index=False, float_format="%.6g")
print(df.to_string(index=False))

# ---- figure ----
fig, ax = plt.subplots(figsize=(10.5, 6.0))
import numpy as np
x = np.arange(len(df))
w = 0.26
labels = [f"{r}\n{g}" for r, g in zip(df["lead_rsid"], df["lead_gene"])]

c_eas, c_eur, c_meta = "#C0392B", "#2C5F8A", "#7F8C8D"  # EAS 红(高频侧), EUR 蓝, meta 灰
b1 = ax.bar(x - w, df["gnomAD_EAS_AF"], w, label="东亚 EAS (gnomAD v4 genome)", color=c_eas)
b2 = ax.bar(x,     df["gnomAD_EUR_NFE_AF"], w, label="欧洲 EUR/NFE (gnomAD v4 genome)", color=c_eur)
b3 = ax.bar(x + w, df["IgAN_meta_EAF"], w, label="跨祖先 meta EAF (GCST90018866)", color=c_meta)

# value labels
for xi, (eas, eur, meta) in enumerate(zip(df["gnomAD_EAS_AF"], df["gnomAD_EUR_NFE_AF"], df["IgAN_meta_EAF"])):
    ax.text(xi - w, eas + 0.015, f"{eas:.3f}", ha="center", fontsize=8.5, color=c_eas)
    ax.text(xi,     eur + 0.015, f"{eur:.3f}", ha="center", fontsize=8.5, color=c_eur)
    if not np.isnan(meta):
        ax.text(xi + w, meta + 0.015, f"{meta:.3f}", ha="center", fontsize=8.5, color=c_meta)

# annotation of Wang2021 P (P-only, direction unknown) & IgAN meta P
for xi, (wang, igan) in enumerate(zip(df["Wang2021_Han_P"], df["IgAN_meta_P"])):
    wang_s = "1.2e-9" if wang < 1e-8 else f"{wang:.3f}"
    note_w = f"Wang2021 汉族 Gd-IgA1 P={wang_s}"
    note_i = f"IgAN meta P={igan:.2f}" if not np.isnan(igan) else "IgAN 文件无 X"
    ax.text(xi, 1.06, note_w, ha="center", fontsize=8.2, color="#8E44AD", weight="bold")
    ax.text(xi, 0.985, note_i, ha="center", fontsize=7.6, color="#555555")

ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=10)
ax.set_ylim(0, 1.12)
ax.set_ylabel("效应/替代等位基因频率 (AF)", fontsize=11)
ax.set_title("M5 东亚异质性：O-糖基化轴 lead 位点的人群频率分层\n(gnomAD v4 genome; EUR 以 NFE 计; rs5910940 为 X 非 PAR 口径)", fontsize=11.5)
ax.legend(loc="upper right", fontsize=9, frameon=False)
ax.spines[["top", "right"]].set_visible(False)
ax.grid(axis="y", ls=":", alpha=0.4)
fig.text(0.01, 0.015, "来源：gnomAD GraphQL API (gnomad_r4/GRCh38) 2026-09-03 实测；跨祖先 EAF 与 P 值见 lead_x_东亚证据汇总.tsv。"
                      "频率为描述性分层信息，正式人群差异检验需等位基因方向一致的 GWAS 汇总数据。", fontsize=7.6, color="#666666")

import os
os.makedirs("figures", exist_ok=True)
fig.savefig("figures/M5_东亚频率分层图.png", dpi=300, bbox_inches="tight")
fig.savefig("figures/M5_东亚频率分层图.pdf", bbox_inches="tight")
print("saved figures/M5_东亚频率分层图.png/.pdf")
