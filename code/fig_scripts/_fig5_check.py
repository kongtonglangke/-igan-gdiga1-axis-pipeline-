"""_fig5_check.py — Fig5 布局自检（text/patch/containment/ticklabels/legend 五件套）

Fig5 = 2 行 1 列上下布局（panel a 上 / panel b 下），legend 为
panel_a/panel_b 函数内 ax.legend（a 右上图内 / b 左下图内）。
GridSpec/figsize 必须与 fig5_render.py main() 一致：
  plt.figure(figsize=(220,205)mm) + GridSpec(2,1,hspace=0.42,height_ratios=[0.92,1.0])
legend 检查：ax.get_legend() bbox vs ax 全内容（text+patch+tick）求交。
"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import sys
sys.path.insert(0, _paths.legacy("阶段4/主图/scripts"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import fig5_render as F
from style_common import mm_to_in, DPI_DRAFT


def iw(r1, r2):
    return min(r1.x1, r2.x1) - max(r1.x0, r2.x0), min(r1.y1, r2.y1) - max(r1.y0, r2.y0)


def check_text(ax, name):
    n = 0; over = 0
    all_t = list(ax.texts) + ([ax.title] if ax.title.get_text() else []) \
        + ([ax.xaxis.label] if ax.xaxis.label.get_text() else []) \
        + ([ax.yaxis.label] if ax.yaxis.label.get_text() else [])
    for i in range(len(all_t)):
        for j in range(i + 1, len(all_t)):
            r1 = all_t[i].get_window_extent()
            r2 = all_t[j].get_window_extent()
            w, h = iw(r1, r2)
            if w <= 1 or h <= 1:
                continue
            over += 1
            a = all_t[i].get_text().replace("\n", " ")[:25]
            b = all_t[j].get_text().replace("\n", " ")[:25]
            print(f"  TXT {w:5.0f}x{h:3.0f}px | \"{a}\" <-> \"{b}\"")
    n = len(all_t)
    print(f"=== panel {name}: {n} texts, {over} text-text overlaps")
    return over


def check_patch(ax, name):
    n = 0; over = 0
    for i in range(len(ax.patches)):
        for j in range(i + 1, len(ax.patches)):
            p1, p2 = ax.patches[i], ax.patches[j]
            # 跳过 facecolor='none' 轮廓框配对（REVERSED 红框/annotate
            # 框与柱 bbox 相交是设计意图，非意外重叠）
            try:
                fc1 = p1.get_facecolor()
                fc2 = p2.get_facecolor()
            except Exception:
                continue
            if (len(fc1) == 4 and fc1[3] == 0) or (len(fc2) == 4 and fc2[3] == 0):
                continue
            r1 = p1.get_window_extent()
            r2 = p2.get_window_extent()
            w, h = iw(r1, r2)
            if w <= 6 or h <= 6:
                continue
            over += 1
            print(f"  PATCH {w:5.0f}x{h:3.0f}px | p#{i} <-> p#{j}")
    n = len(ax.patches)
    print(f"=== panel {name}: {n} patches, {over} patch-patch overlaps")
    return over


def check_containment(ax, name):
    leaks = 0
    texts = list(ax.texts)
    for t in texts:
        r = t.get_window_extent()
        cx, cy = (r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2
        for p in ax.patches:
            pr = p.get_window_extent()
            if pr.x0 <= cx <= pr.x1 and pr.y0 <= cy <= pr.y1:
                if not (pr.x0 <= r.x0 and r.x1 <= pr.x1
                        and pr.y0 <= r.y0 and r.y1 <= pr.y1):
                    if (r.x0 < pr.x0 - 3 or r.x1 > pr.x1 + 3
                            or r.y0 < pr.y0 - 3 or r.y1 > pr.y1 + 3):
                        leaks += 1
                        print(f"  ESCAPE {r.width:.0f}x{r.height:.0f}px | "
                              f"\"{t.get_text()[:25]}\" overflow box")
    print(f"=== panel {name}: containment leaks={leaks}")
    return leaks


def check_ticklabels(ax, name):
    n = 0; over = 0
    all_ticks = list(ax.get_xticklabels()) + list(ax.get_yticklabels())
    TICK_TOL = 2
    for i in range(len(all_ticks)):
        for j in range(i + 1, len(all_ticks)):
            r1 = all_ticks[i].get_window_extent()
            r2 = all_ticks[j].get_window_extent()
            w, h = iw(r1, r2)
            if w <= TICK_TOL or h <= TICK_TOL:
                continue
            over += 1
            a = all_ticks[i].get_text().replace("\n", " ")[:20]
            b = all_ticks[j].get_text().replace("\n", " ")[:20]
            print(f"  XTICK {w:5.0f}x{h:3.0f}px | \"{a}\" <-> \"{b}\"")
    n = len(all_ticks)
    print(f"=== panel {name}: {n} ticklabels, {over} tick-tick overlaps")
    return over


def check_legend_overlap(ax, name, tol=3):
    """ax 内 legend bbox vs ax 全内容（text+patch+tick）求交。
    legend 自身文字/框不在 ax.texts/ax.patches 内，不会自交。"""
    over = 0
    leg = ax.get_legend()
    if leg is None:
        print(f"=== panel {name}: no legend, legend-overlap=0")
        return 0
    lr = leg.get_window_extent()
    items = (list(ax.texts) + list(ax.patches)
             + list(ax.get_xticklabels()) + list(ax.get_yticklabels())
             + ([ax.title] if ax.title.get_text() else [])
             + ([ax.xaxis.label] if ax.xaxis.label.get_text() else [])
             + ([ax.yaxis.label] if ax.yaxis.label.get_text() else []))
    for it in items:
        r = it.get_window_extent()
        w, h = iw(lr, r)
        if w > tol and h > tol:
            over += 1
            s = it.get_text().replace("\n", " ")[:20] if hasattr(it, "get_text") else str(it)
            print(f"  LEG-OVER {w:.0f}x{h:.0f}px | legend <-> \"{s}\"")
    print(f"=== panel {name}: legend-overlap={over}")
    return over


def main():
    fig = plt.figure(figsize=(mm_to_in(220), mm_to_in(205)), dpi=DPI_DRAFT)
    gs = GridSpec(2, 1, figure=fig, height_ratios=[0.92, 1.0], hspace=0.42)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[1, 0])
    F.panel_a(ax_a)
    F.panel_b(ax_b)
    fig.canvas.draw()

    total = 0
    for k, ax in enumerate([ax_a, ax_b]):
        print(f"\n--- panel {chr(97 + k)} ---")
        total += check_text(ax, chr(97 + k))
        total += check_patch(ax, chr(97 + k))
        total += check_containment(ax, chr(97 + k))
        total += check_ticklabels(ax, chr(97 + k))
        total += check_legend_overlap(ax, chr(97 + k))
    print(f"\nTOTAL overlaps: {total}")


if __name__ == "__main__":
    main()
