"""style_common.py — 阶段4 主图统一视觉规范（GM 口径 + Okabe-Ito 色盲安全板）
作者：AI 王严课题组代办｜2026-09-06｜依赖：matplotlib ≥3.8, pandas ≥2, numpy ≥2

引用：fig43 venv（阶段4 绘图专用隔离环境）。
"""
from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# === 期刊硬约束（Springer/Genome Medicine）===
DPI_DRAFT = 150           # v01 预览 dpi（便于快速查看）
DPI_TIFF = 600            # 投稿版线/组合 dpi
FIG_WIDTH_MM = 183        # GM 双栏宽度（mm）；单栏为 89 mm
FIG_HEIGHT_MM_BY_PANEL = 55  # 每面板推荐高度（mm），按内容可调

# === 字体（按可用性回退；Arial 未装时退 DejaVu Sans）===
def _setup_fonts():
    # 优先 DejaVu Sans（Linux/Win 通用，Unicode 全覆盖，包括 ⁻⁸⁶² 等上标负号）；
    # 其次 Arial（投稿正文首选无衬线），再次 Helvetica。
    candidates = ["DejaVu Sans", "Arial", "Helvetica", "Liberation Sans"]
    chosen = "DejaVu Sans"
    for name in candidates:
        try:
            mpl.font_manager.findfont(name, fallback_to_default=False)
            chosen = name
            break
        except Exception:  # noqa: BLE001
            continue
    mpl.rcParams["font.family"] = ["sans-serif"]
    mpl.rcParams["font.sans-serif"] = [chosen] + [c for c in candidates if c != chosen]
    mpl.rcParams["pdf.fonttype"] = 42   # 嵌入 TrueType（投稿 PDF 必备）
    mpl.rcParams["ps.fonttype"] = 42
    mpl.rcParams["axes.unicode_minus"] = False
    mpl.rcParams["mathtext.fontset"] = "dejavusans"

_setup_fonts()

# === 字号（双栏缩放后的等效 pt；实际字号将按宽度再缩放）===
FS_PANEL_LABEL = 12       # 面板字母（a/b/c…），粗体
FS_AXIS_LABEL = 9.5       # 轴标
FS_TICK = 8.5             # 刻度
FS_LEGEND = 8.5           # 图例
FS_INCELL = 7.5           # 矩阵格内数字
FS_TITLE = 11             # 面板标题（如有）

# === Okabe-Ito 色盲安全 8 色 + 中性灰 ===
OKABE_ITO = {
    "black":      "#000000",
    "orange":     "#E69F00",
    "sky_blue":   "#56B4E9",
    "bluish_grn": "#009E73",
    "yellow":     "#F0E442",
    "blue":       "#0072B2",
    "vermillion": "#D55E00",
    "reddish_pur": "#CC79A7",
    "gray":       "#999999",
}
# 课题约定两组对比：Kiryluk'17 rs13226913 ↔ Wang'21 rs10238682
C_RS13226913 = OKABE_ITO["blue"]        # 蓝
C_RS10238682 = OKABE_ITO["vermillion"]  # 朱红
C_RS7856182 = OKABE_ITO["orange"]       # 橙
C_RS5910940 = OKABE_ITO["bluish_grn"]   # 蓝绿
C_NOTTEST = OKABE_ITO["gray"]
C_GW_LINE = OKABE_ITO["vermillion"]     # genome-wide sig 红虚线

# 热图主色阶（sequential 单色，从白→深红，色盲安全）
HEAT_CMAP = mpl.colors.LinearSegmentedColormap.from_list(
    "okabe_heat_white_red",
    ["#FFFFFF", "#FDE0DD", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"],
)

# === 灰格 + 斜纹（"not testable"）===
def hatch_not_tested(ax, x, y, label="NA", text_color=None):
    """在 (x,y) 处画灰底斜纹方块并加 NA 文字。"""
    ax.add_patch(plt.Rectangle((x - 0.5, y - 0.5), 1, 1,
                               facecolor="#F0F0F0",
                               edgecolor=OKABE_ITO["gray"],
                               hatch="///",
                               lw=0.8, zorder=0))
    ax.text(x, y, label, ha="center", va="center",
            fontsize=FS_INCELL - 0.5,
            color=text_color or OKABE_ITO["gray"], zorder=3)

# === 高亮方框（红框）===
def highlight_box(ax, x0, x1, y0, y1, color=None, lw=1.8, ls="-"):
    color = color or C_GW_LINE
    ax.plot([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0],
            color=color, lw=lw, ls=ls, zorder=5)

# === 工具：mm → inches ===
def mm_to_in(mm: float) -> float:
    return mm / 25.4

# === 工具：创建 3 面板垂直堆叠 Figure（双栏宽）===
def new_fig_three_panels(panel_h_mm=(60, 70, 65), w_mm=FIG_WIDTH_MM, dpi=DPI_DRAFT):
    """返回 (fig, axes)；axes 长度 3，垂直堆叠。使用 constrained_layout 自动间距。"""
    total_h_mm = sum(panel_h_mm) + 14
    fig, axes = plt.subplots(
        3, 1,
        figsize=(mm_to_in(w_mm), mm_to_in(total_h_mm)),
        dpi=dpi,
        constrained_layout=True,
        gridspec_kw={"height_ratios": panel_h_mm, "hspace": 0.35},
    )
    return fig, axes

# === 公共保存 ===
SAVE_BBOX_KW = dict(bbox_inches="tight", pad_inches=0.4)

def save_v01(fig, out_path: str):
    fig.savefig(out_path, dpi=DPI_DRAFT, **SAVE_BBOX_KW)
    fig.savefig(out_path.replace(".png", ".pdf"), dpi=DPI_DRAFT, **SAVE_BBOX_KW)


def save_final(fig, out_stem: str, dpi_tiff=DPI_TIFF):
    """投稿版导出：PNG 300 dpi（预览/审稿）+ PDF（矢量）+ TIFF 600 dpi（Springer 投稿首选）。

    out_stem 不含扩展名（如 .../Fig2_v1.0）。TIFF 用 LZW 压缩 + 白底（组合图无透明）。
    所有格式统一 bbox_inches='tight' + pad_inches=0.4：避免边缘文字/坐标轴被裁切，
    文字/轴数字不重叠画框。
    """
    fig.savefig(out_stem + ".png", dpi=300, **SAVE_BBOX_KW)
    fig.savefig(out_stem + ".pdf", dpi=300, **SAVE_BBOX_KW)
    try:
        fig.savefig(out_stem + ".tiff", dpi=dpi_tiff,
                    pil_kwargs={"compression": "tiff_lzw"}, facecolor="white",
                    **SAVE_BBOX_KW)
    except Exception as e:  # noqa: BLE001
        print(f"TIFF save failed ({e}); retrying without compression")
        fig.savefig(out_stem + ".tiff", dpi=dpi_tiff, facecolor="white",
                    **SAVE_BBOX_KW)