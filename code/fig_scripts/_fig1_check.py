# -*- coding: utf-8 -*-
"""_fig1_check.py — Fig1 布局自检：渲染后测量 text/patch 像素边界框，
输出重叠对与每行位置，作为"无重叠"验收的客观依据（替代人工目检）。"""
import os, sys
os.environ['OUT_STEM'] = 'Fig1_chk'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fig1_render as F
from matplotlib.gridspec import GridSpec

fig = plt.figure(figsize=(F.mm_to_in(183), F.mm_to_in(280)), dpi=F.DPI_DRAFT)
gs = GridSpec(3, 1, figure=fig, height_ratios=[100, 62, 108], hspace=0.42)
axs = [fig.add_subplot(gs[i, 0]) for i in range(3)]
F.panel_a(axs[0]); F.panel_b(axs[1]); F.panel_c(axs[2])
fig.canvas.draw()


def skip_empty(t):
    return t.get_text().strip() == ""


def bbox_area(r):
    return max(0.0, r.width * r.height)


def inter_wh(r1, r2):
    x0 = max(r1.x0, r2.x0)
    x1 = min(r1.x1, r2.x1)
    y0 = max(r1.y0, r2.y0)
    y1 = min(r1.y1, r2.y1)
    return max(0.0, x1 - x0), max(0.0, y1 - y0)


def check_text(ax, name):
    texts = [t for t in ax.texts if not skip_empty(t)]
    n = len(texts)
    over = 0
    for i in range(n):
        ri = texts[i].get_window_extent()
        for j in range(i + 1, n):
            rj = texts[j].get_window_extent()
            w, h = inter_wh(ri, rj)
            if w <= 1 or h <= 1:
                continue
            area = w * h
            over += 1
            a = texts[i].get_text().replace("\n", " ")[:40]
            b = texts[j].get_text().replace("\n", " ")[:40]
            print(f"  TXT  {area:8.0f}px2 | \"{a}\" <-> \"{b}\"")
    print(f"=== panel {name}: {n} texts, {over} text-text overlaps")
    return over


def check_patch(ax, name):
    """检测 patch 两两重叠（除同色背景填充带外）。"""
    from matplotlib.patches import Rectangle, FancyBboxPatch
    ps = [p for p in ax.patches if isinstance(p, (Rectangle, FancyBboxPatch))]
    n = len(ps)
    over = 0
    for i in range(n):
        ri = ps[i].get_window_extent()
        for j in range(i + 1, n):
            rj = ps[j].get_window_extent()
            w, h = inter_wh(ri, rj)
            if w <= 6 or h <= 6:
                continue
            area = w * h
            over += 1
            print(f"  PATCH {area:8.0f}px2 | patch#{i} <-> patch#{j}")
    print(f"=== panel {name}: {n} patches, {over} patch-patch overlaps")
    return over


def check_containment(ax, name):
    """检测文字是否溢出所属框（含 0~3px 容差以容忍 round 圆角与线宽）。
    算法：对每个 text 找其中心点所在的 patch（若有），再检查 text 的 4 边
    是否都在该 patch 内（容差 tol_px）。
    """
    from matplotlib.patches import Rectangle, FancyBboxPatch
    tol_px = 3
    texts = [t for t in ax.texts if not skip_empty(t)]
    ps = [p for p in ax.patches if isinstance(p, (Rectangle, FancyBboxPatch))]
    over = 0
    for ti, t in enumerate(texts):
        if not t.get_text().strip():
            continue
        ctr = t.get_window_extent()  # text bbox
        cx = (ctr.x0 + ctr.x1) / 2
        cy = (ctr.y0 + ctr.y1) / 2
        # 找包含中心的 patch
        host = None
        for pi, p in enumerate(ps):
            pr = p.get_window_extent()
            if pr.x0 - 2 <= cx <= pr.x1 + 2 and pr.y0 - 2 <= cy <= pr.y1 + 2:
                host = (pi, pr)
                break
        if host is None:
            continue
        pi, pr = host
        # 检查 text 4 边是否在 patch 内（含 tol_px）
        outside_l = max(0, pr.x0 + tol_px - ctr.x0)
        outside_r = max(0, ctr.x1 - (pr.x1 - tol_px))
        outside_b = max(0, pr.y0 + tol_px - ctr.y0)
        outside_t = max(0, ctr.y1 - (pr.y1 - tol_px))
        leaks = sum(int(x > 0) for x in (outside_l, outside_r, outside_b, outside_t))
        if leaks == 0:
            continue
        over += leaks
        sides = []
        if outside_l > 0: sides.append(f"left={outside_l:.0f}px")
        if outside_r > 0: sides.append(f"right={outside_r:.0f}px")
        if outside_b > 0: sides.append(f"bottom={outside_b:.0f}px")
        if outside_t > 0: sides.append(f"top={outside_t:.0f}px")
        body = t.get_text().replace("\n", " ")[:36]
        print(f"  ESCAPE  patch#{pi} | \"{body}\" | {' '.join(sides)}")
    print(f"=== panel {name}: containment leaks={over}")
    return over


total = 0
for k, ax in enumerate(axs):
    total += check_text(ax, chr(97 + k))
    total += check_patch(ax, chr(97 + k))
    total += check_containment(ax, chr(97 + k))

# panel c 行位置清单
axc = axs[2]
print("\n--- panel c row text y positions (data coords) ---")
rows_c = axc.texts
for t in rows_c:
    xy = t.get_position()
    txt = t.get_text().replace("\n", " ")[:36]
    print(f"  y={xy[1]:6.2f}  x={xy[0]:5.2f}  \"{txt}\"")

print(f"\nTOTAL overlaps: {total}")
fig.savefig(os.path.join(F.OUT_DIR, "Fig1_chk.png"), dpi=150, bbox_inches="tight", pad_inches=0.4)
print("saved Fig1_chk.png")
