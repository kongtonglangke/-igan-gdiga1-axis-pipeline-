# -*- coding: utf-8 -*-
"""修复版 Figure_S3：复用 stepC_plot.py 图1主体，仅渲染 lead×语境矩阵热条"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import os, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import matplotlib.patches as mpatches
for _f in ("Microsoft YaHei", "SimHei", "SimSun"):
    try:
        font_manager.findfont(_f, fallback_to_default=False)
        plt.rcParams["font.sans-serif"] = [_f, "DejaVu Sans"]; break
    except Exception: continue
plt.rcParams["axes.unicode_minus"] = False

HERE = _paths.at("阶段2/主证据1_lead_eQTL")
DST = _paths.at("阶段4.5_复现包/submission/additional_file_1")
LEADS = ["rs13226913", "rs10238682", "rs7856182", "rs5910940"]
LEAD_LABEL = {"rs13226913": "rs13226913\n(C1GALT1, Kiryluk'17)",
              "rs10238682": "rs10238682\n(C1GALT1, Wang'21)",
              "rs7856182": "rs7856182\n(GALNT12, +19 kb, Wang'21)",
              "rs5910940": "rs5910940\n(C1GALT1C1 intronic, Kiryluk'17)"}
DSET_ORDER = ["OneK1K_B_naive", "OneK1K_B_memory", "OneK1K_B_intermediate",
              "GTEx_v10_blood", "GTEx_v10_kidney_cortex"]
DSET_LABEL = {"OneK1K_B_naive": "OneK1K B naive\n(n=876)",
              "OneK1K_B_memory": "OneK1K B memory\n(n=726)",
              "OneK1K_B_intermediate": "OneK1K B inter.\n(n=618)",
              "GTEx_v10_blood": "GTEx v10 whole blood\n(n=853)",
              "GTEx_v10_kidney_cortex": "GTEx v10 kidney cortex\n(n=106)"}

rows = [l.rstrip("\n").split("\t") for l in open(HERE("lead_eQTL_matrix.tsv"), encoding="utf-8")]
rows = [r for r in rows if r and not r[0].startswith("lead_rsid")]
def get(lead, ds):
    for r in rows:
        if r[0] == lead and r[5] == ds: return r
    return None

fig, ax = plt.subplots(figsize=(9.5, 4.6))
X = np.arange(len(DSET_ORDER)); width = 0.19
colors_ds = {"OneK1K_B_naive": "#2166ac", "OneK1K_B_memory": "#4393c3",
             "OneK1K_B_intermediate": "#92c5de", "GTEx_v10_blood": "#d6604d",
             "GTEx_v10_kidney_cortex": "#b2182b"}
handles = []
for li, lead in enumerate(LEADS):
    xpos = X + (li - 1.5) * width
    for xi, ds in enumerate(DSET_ORDER):
        r = get(lead, ds); axv = xpos[xi]
        if r is None: continue
        ap = r[9]
        if ap == "LEAD_ABSENT":
            ax.bar(axv, 1.0, width, color="#969696", edgecolor="k", lw=0.4)
            # 不在每个柱顶标字，避免重叠；用下方 group 计数统一说
        elif ap == "":
            ax.bar(axv, 1.0, width, color="#d9d9d9", edgecolor="k", lw=0.4, hatch="//")
            # 不在每个柱顶标字，避免重叠；用下方 group 计数统一说
        else:
            p = float(ap); lp = min(-np.log10(p), 12)
            ax.bar(axv, lp, width, color=colors_ds[ds], edgecolor="k", lw=0.4)
            sig = "***" if p < 5e-8 else ("**" if p < 1e-4 else ("*" if p < 0.05 else ""))
            if sig:
                ax.text(axv, lp + 0.25, sig, ha="center", va="bottom", fontsize=9, fontweight="bold")
for ds in DSET_ORDER:
    handles.append(mpatches.Patch(color=colors_ds[ds], label=DSET_LABEL[ds]))
handles += [mpatches.Patch(color="#969696", label="lead SNP absent from dataset"),
            mpatches.Patch(color="#d9d9d9", hatch="//", label="annotated gene not detected (low B-cell expression)")]
ax.axhline(-np.log10(5e-8), color="red", ls="--", lw=0.9)
ax.text(-0.20, -np.log10(5e-8) + 0.25, "genome-wide sig (5e-8)", color="red", fontsize=7.5, ha="left")
ax.axhline(-np.log10(1e-4), color="gray", ls=":", lw=0.8)
ax.set_xticks(X)
# 在每个 X context 下方显示 n_SNP_absent / n_gene_not_tested 计数（防柱顶标签重叠）
_per_x_stats = {xi: {"absent": 0, "untested": 0} for xi in range(len(DSET_ORDER))}
for li, lead in enumerate(LEADS):
    for xi, ds in enumerate(DSET_ORDER):
        r = get(lead, ds)
        if r is None: continue
        _ap = r[9]
        if _ap == "LEAD_ABSENT": _per_x_stats[xi]["absent"] += 1
        elif _ap == "": _per_x_stats[xi]["untested"] += 1
for xi, st in _per_x_stats.items():
    parts = []
    if st["absent"]: parts.append(f"{st['absent']}\u00d7SNP absent")
    if st["untested"]: parts.append(f"{st['untested']}\u00d7gene not tested")
    if parts:
        # v3.3：y 从 -2.4 → -2.1（ax 内）；X 轴 tick "n=xxx" 在 y≈-1.0，
        # 让统计文字下移到 -2.1 仍与 x tick 有 1.1 unit (~30px) 间距
        ax.text(xi, -2.1, ", ".join(parts), ha="center", va="top",
                fontsize=6.5, color="#444444", style="italic")
ax.set_xticklabels([DSET_LABEL[d] for d in DSET_ORDER], fontsize=8.5)
ax.tick_params(axis="x", pad=4)
ax.set_ylabel("−log10(p) (minimal nominal P of lead SNP on the annotated gene)", fontsize=9.5)
ax.set_ylim(0, 14.5)
ax.set_title("cis-eQTL attribution of the four Gd-IgA1 lead SNPs across contexts\n(B-cell / plasma-cell mechanism)", fontsize=11.5)
# 图例从 ax 内右上角挪到画布外下方，避免遮挡 GTEx v10 kidney cortex 高柱与星号
fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.005),
           ncol=3, fontsize=7.5, frameon=True, framealpha=0.92, edgecolor="#888888")
fig.subplots_adjust(left=0.08, right=0.97, top=0.90, bottom=0.34)
stem = DST("Figure_S3_lead_eQTL_full_matrix_heat")
fig.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
try:
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
except PermissionError as e:
    print("[warn] PDF 写入被外部进程占用（Permission denied），PNG 已写入: {stem}.pdf 待用户关闭占用进程后重渲染 PDF", e)
plt.close(fig)
print(f"saved {stem}.png/.pdf")