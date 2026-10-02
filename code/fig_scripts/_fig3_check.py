"""_fig3_check.py — Fig3 布局自检（text/patch/containment/ticklabels 四件套）

输出每个面板的客观重叠清单（像素级窗口 bbox 两两求交）。
GridSpec/figsize 必须与 fig3_render.py 一致（new_fig_three_panels(panel_h_mm=(60,76,60))）。
"""
from __future__ import annotations

import os
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import fig3_render as F


def iw(r1, r2):
    x0 = max(r1.x0, r2.x0); x1 = min(r1.x1, r2.x1)
    y0 = max(r1.y0, r2.y0); y1 = min(r1.y1, r2.y1)
    return max(0.0, x1 - x0), max(0.0, y1 - y0)


def skip_empty(t):
    return not (t.get_text() or "").strip()


def check_text(ax, name):
    texts = [t for t in ax.texts if not skip_empty(t)]
    n = len(texts)
    over = 0
    for i in range(n):
        ri = texts[i].get_window_extent()
        for j in range(i + 1, n):
            rj = texts[j].get_window_extent()
            w, h = iw(ri, rj)
            if w <= 1 or h <= 1:
                continue
            over += 1
            a = texts[i].get_text().replace("\n", " ")[:40]
            b = texts[j].get_text().replace("\n", " ")[:40]
            print(f"  TXT  {w:5.0f}x{h:3.0f}px | \"{a}\" <-> \"{b}\"")
    print(f"=== panel {name}: {n} texts, {over} text-text overlaps")
    return over


def check_patch(ax, name):
    from matplotlib.patches import Rectangle, FancyBboxPatch
    ps = [p for p in ax.patches if isinstance(p, (Rectangle, FancyBboxPatch))]
    n = len(ps)
    over = 0
    for i in range(n):
        ri = ps[i].get_window_extent()
        for j in range(i + 1, n):
            rj = ps[j].get_window_extent()
            w, h = iw(ri, rj)
            if w <= 6 or h <= 6:
                continue
            over += 1
            print(f"  PATCH {w:5.0f}x{h:3.0f}px | patch#{i} <-> patch#{j}")
    print(f"=== panel {name}: {n} patches, {over} patch-patch overlaps")
    return over


def check_containment(ax, name):
    from matplotlib.patches import Rectangle, FancyBboxPatch
    # tol_px=4 容忍 matplotlib font metric descender padding (~3 px) +
    # lw=0.3 edgecolor 边沿 (0.5 px)，与视觉一致；之前 3px 误报 PPH 堆叠
    # 条内 n=xxxx 文字 2px 边沿 leak。
    tol_px = 4
    texts = [t for t in ax.texts if not skip_empty(t)]
    ps = [p for p in ax.patches if isinstance(p, (Rectangle, FancyBboxPatch))]
    over = 0
    for t in texts:
        ctr = t.get_window_extent()
        cx = (ctr.x0 + ctr.x1) / 2; cy = (ctr.y0 + ctr.y1) / 2
        host = None
        for pi, p in enumerate(ps):
            pr = p.get_window_extent()
            if pr.x0 - 2 <= cx <= pr.x1 + 2 and pr.y0 - 2 <= cy <= pr.y1 + 2:
                host = (pi, pr); break
        if host is None:
            continue
        pi, pr = host
        ol = max(0, pr.x0 + tol_px - ctr.x0)
        orr = max(0, ctr.x1 - (pr.x1 - tol_px))
        ob = max(0, pr.y0 + tol_px - ctr.y0)
        ot = max(0, ctr.y1 - (pr.y1 - tol_px))
        # 用 tol_px 阈值：matplotlib font descender padding 让 text bbox
        # 比可见字符大 ~3 px，与视觉一致的小 leak 不计入 leak 数量。
        leaks = sum(int(x > tol_px) for x in (ol, orr, ob, ot))
        if leaks == 0:
            continue
        over += leaks
        sides = []
        if ol > 0: sides.append(f"left={ol:.0f}px")
        if orr > 0: sides.append(f"right={orr:.0f}px")
        if ob > 0: sides.append(f"bottom={ob:.0f}px")
        if ot > 0: sides.append(f"top={ot:.0f}px")
        print(f"  ESCAPE patch#{pi} | \"{t.get_text().replace(chr(10), ' ')[:36]}\" | {' '.join(sides)}")
    print(f"=== panel {name}: containment leaks={over}")
    return over


def check_ticklabels(ax, name):
    n = 0; over = 0
    all_ticks = list(ax.get_xticklabels()) + list(ax.get_yticklabels())
    # TICK_TOL=3：matplotlib font metric 让 y tick bbox y 高度比可见字
    # 高 ~3 px；跨类型 tick（x-tick 在 ax 底 vs y-tick 在 ax 左）bbox
    # 物理分离但跨类型比较时 y 方向边沿仍可能 2.4~2.9 px 微触，与视
    # 觉无感对齐。TICK_TOL=2 时 panel b 报 9 个跨类型假阳性；放宽到 3
    # 让其跳过。Panel c 真实 tick-tick 重叠（如 v3.2-rotate 42×21px）
    # h=21>>3 仍报。
    TICK_TOL = 3
    for i in range(len(all_ticks)):
        for j in range(i + 1, len(all_ticks)):
            r1 = all_ticks[i].get_window_extent()
            r2 = all_ticks[j].get_window_extent()
            w, h = iw(r1, r2)
            if w <= TICK_TOL or h <= TICK_TOL:
                continue
            over += 1
            a = all_ticks[i].get_text().replace("\n", " ")[:30]
            b = all_ticks[j].get_text().replace("\n", " ")[:30]
            print(f"  XTICK {w:5.0f}x{h:3.0f}px | \"{a}\" <-> \"{b}\"")
    n = len(all_ticks)
    print(f"=== panel {name}: {n} ticklabels, {over} tick-tick overlaps")
    return over


def check_legend_overlap(ax, name):
    """检查 ax.get_legend() 框与 ax.patches 的窗口重叠。

    matplotlib legend 是独立 Artist 容器（ax.legend_），其内部 handle 与
    frame 的 bbox 不在 ax.patches / ax.texts 列表中，原四件套检测不到
    "legend 框盖住柱子或文字"。本函数枚举 ax.patches 的 window extent 与
    legend window extent 求交，TOL=2px 容忍 legend edge 1px + 安全边沿。
    若 ax.get_legend() 返回 None（无 legend）则跳过。
    """
    from matplotlib.patches import Rectangle, FancyBboxPatch
    leg = ax.get_legend()
    if leg is None:
        print(f"=== panel {name}: no legend, skipped")
        return 0
    rl = leg.get_window_extent()
    ps = [p for p in ax.patches if isinstance(p, (Rectangle, FancyBboxPatch))]
    over = 0
    # 关联 handles（每行一个 handle patch）
    handles = [h for h in leg.legend_handles if h is not None] if hasattr(leg, "legend_handles") else []
    for i, p in enumerate(ps):
        rp = p.get_window_extent()
        # legend 整体框 vs patch
        w, h = iw(rl, rp)
        if w > 2 and h > 2:
            over += 1
            print(f"  LEG-FRAME {w:5.0f}x{h:3.0f}px | legend<->patch#{i} (legend box covers patch)")
        # legend 关联 handles vs patch（handles patch 在 legend 框内对齐，故 vs 数据 patch
        # 重叠等价于 legend-frame vs patch 重叠；此处略避免重复计数）
    print(f"=== panel {name}: legend overlap with {len(ps)} patches = {over}")
    return over


fig, axs = F.new_fig_three_panels(panel_h_mm=(60, 76, 60))
F.panel_a(axs[0]); F.panel_b(axs[1]); F.panel_c(axs[2])
fig.canvas.draw()

total = 0
for k, ax in enumerate(axs):
    print(f"\n--- panel {chr(97 + k)} ---")
    total += check_text(ax, chr(97 + k))
    total += check_patch(ax, chr(97 + k))
    total += check_containment(ax, chr(97 + k))
    total += check_ticklabels(ax, chr(97 + k))
    total += check_legend_overlap(ax, chr(97 + k))
print(f"\nTOTAL overlaps: {total}")