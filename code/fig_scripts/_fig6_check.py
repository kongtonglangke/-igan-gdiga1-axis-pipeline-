"""_fig6_check.py — Fig6 布局自检（text/patch/containment/ticklabels/legend 五件套）

Fig6 = 2×2 网格（a 堆叠条 / b 热图 / c 散点 / d 示意）+ panel a 的
fig-level legend（fig.legends）。panel b 另有 ax.legend（右侧外置）。
GridSpec/figsize 必须与 fig6_render.py main() 一致：
  plt.figure(figsize=(220,180)mm)
  GridSpec(2,2,height_ratios=[1,1.20],width_ratios=[1,1],hspace=0.50,wspace=0.40)
legend 检查：fig.legends × 每 ax 全内容 + ax.get_legend() × 自身 ax 内容。
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
import fig6_render as F
from style_common import mm_to_in, DPI_DRAFT


def iw(r1, r2):
    return min(r1.x1, r2.x1) - max(r1.x0, r2.x0), min(r1.y1, r2.y1) - max(r1.y0, r2.y0)


def check_text(ax, name):
    n = 0; over = 0
    all_t = []
    for t in list(ax.texts):
        if t.get_text().strip():          # 跳过 annotate 空文本（箭头伪 text）
            all_t.append(t)
    all_t += ([ax.title] if ax.title.get_text() else []) \
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
            # 跳过 facecolor='none' / alpha=0 轮廓框配对（设计意图：
            # Fig5 REVERSED 红框、panel d 边框等与柱 bbox 相交）
            try:
                fc1 = p1.get_facecolor()
                fc2 = p2.get_facecolor()
            except Exception:
                continue
            if (len(fc1) == 4 and fc1[3] == 0) or (len(fc2) == 4 and fc2[3] == 0):
                continue
            # 跳过 #FFCCCC 风险扣分斜纹块配对（panel a：-P 斜纹覆盖 F 段
            # 尾部是"扣除"语义设计意图，非意外重叠；两处 44×34px）
            def _is_penalty(fc):
                return (len(fc) == 4 and round(fc[0], 2) == 1.00
                        and round(fc[1], 2) == 0.80 and round(fc[2], 2) == 0.80)
            if _is_penalty(fc1) or _is_penalty(fc2):
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
    texts = [t for t in ax.texts if t.get_text().strip()]
    for t in texts:
        r = t.get_window_extent()
        cx, cy = (r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2
        for p in ax.patches:
            pr = p.get_window_extent()
            # 跳过 axvspan（alpha<1 的半透明背景 patch）：axvspan 默认
            # 占据整个 ax ylim 高度，所有 text 中心 y 都在其内，bbox 越界
            # 是 text 宽度超出 axvspan x 范围，与"文字溢出封闭 box"语义不符。
            try:
                fc = p.get_facecolor()
                if len(fc) == 4 and fc[3] < 0.95:
                    continue
            except Exception:
                pass
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


def _ax_items(ax):
    return (list(ax.texts) + list(ax.patches)
            + list(ax.get_xticklabels()) + list(ax.get_yticklabels())
            + ([ax.title] if ax.title.get_text() else [])
            + ([ax.xaxis.label] if ax.xaxis.label.get_text() else [])
            + ([ax.yaxis.label] if ax.yaxis.label.get_text() else []),
            [ax.title] if ax.title.get_text() else [])


def check_legend_overlap(fig, axs, tol=3):
    """fig.legends + 每 ax 的 ax.get_legend()，bbox vs 全 ax 内容求交。"""
    over = 0
    legends = list(fig.legends)
    for ax in axs:
        leg = ax.get_legend()
        if leg is not None:
            legends.append(leg)
    for leg in legends:
        lr = leg.get_window_extent()
        lname = "FIG-LEG" if leg in fig.legends else "AX-LEG"
        for k, ax in enumerate(axs):
            items, _ = _ax_items(ax)
            for it in items:
                r = it.get_window_extent()
                w, h = iw(lr, r)
                if w > tol and h > tol:
                    over += 1
                    s = it.get_text().replace("\n", " ")[:20] if hasattr(it, "get_text") else str(it)
                    print(f"  LEG-OVER {w:.0f}x{h:.0f}px | {lname} <-> panel {chr(97 + k)} \"{s}\"")
    print(f"=== legends ({len(legends)}): legend-overlap={over}")
    return over


def main():
    fig = plt.figure(figsize=(mm_to_in(220), mm_to_in(180)), dpi=DPI_DRAFT)
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1, 1.20],
                  width_ratios=[1, 1], hspace=0.50, wspace=0.40)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])
    F.panel_a(ax_a)
    F.panel_b(ax_b)
    F.panel_c(ax_c)
    F.panel_d(ax_d)
    # v3.1：移除 fig-level legend，panel a/b 的 ax.legend 在各自 panel 内
    # (ax_a 内右下 + ax_b 外右上)。fig-level 仅检测空集合。
    fig.canvas.draw()

    axs = [ax_a, ax_b, ax_c, ax_d]
    total = 0
    for k, ax in enumerate(axs):
        print(f"\n--- panel {chr(97 + k)} ---")
        total += check_text(ax, chr(97 + k))
        total += check_patch(ax, chr(97 + k))
        total += check_containment(ax, chr(97 + k))
        total += check_ticklabels(ax, chr(97 + k))
    total += check_legend_overlap(fig, axs)
    print(f"\nTOTAL overlaps: {total}")


if __name__ == "__main__":
    main()
