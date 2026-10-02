# -*- coding: utf-8 -*-
"""_fig2_check.py — Fig2 布局自检 + panel b 柱颜色像素核对。
自检部分：check_text/check_patch/check_containment（复用 fig1 骨架）。
颜色核对：canvas buffer 采样，panel b 每根柱中心 RGB vs 期望 lead 颜色，
  捕获"第一个柱子颜色混了"类问题（柱间距 0 → 间隙、色块贴边/重叠）。
"""
import os, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fig2_render as F

fig, axs = F.new_fig_three_panels(panel_h_mm=(60, 70, 60))
F.panel_a(axs[0]); F.panel_b(axs[1]); F.panel_c(axs[2])
fig.canvas.draw()


def skip_empty(t):
    return t.get_text().strip() == ""


def inter_wh(r1, r2):
    x0 = max(r1.x0, r2.x0); x1 = min(r1.x1, r2.x1)
    y0 = max(r1.y0, r2.y0); y1 = min(r1.y1, r2.y1)
    return max(0.0, x1 - x0), max(0.0, y1 - y0)


def check_text(ax, name):
    texts = [t for t in ax.texts if not skip_empty(t)]
    n = len(texts); over = 0
    for i in range(n):
        ri = texts[i].get_window_extent()
        for j in range(i + 1, n):
            rj = texts[j].get_window_extent()
            w, h = inter_wh(ri, rj)
            if w <= 1 or h <= 1:
                continue
            over += 1
            a = texts[i].get_text().replace("\n", " ")[:36]
            b = texts[j].get_text().replace("\n", " ")[:36]
            print(f"  TXT {w*h:8.0f}px2 | \"{a}\" <-> \"{b}\"")
    print(f"=== panel {name}: {n} texts, {over} text-text overlaps")
    return over


def check_patch(ax, name):
    from matplotlib.patches import Rectangle, FancyBboxPatch
    ps = [p for p in ax.patches if isinstance(p, (Rectangle, FancyBboxPatch))]
    n = len(ps); over = 0
    for i in range(n):
        ri = ps[i].get_window_extent()
        for j in range(i + 1, n):
            rj = ps[j].get_window_extent()
            w, h = inter_wh(ri, rj)
            if w <= 6 or h <= 6:
                continue
            over += 1
            print(f"  PATCH {w*h:8.0f}px2 | patch#{i} <-> patch#{j}")
    print(f"=== panel {name}: {n} patches, {over} patch-patch overlaps")
    return over


def check_containment(ax, name):
    from matplotlib.patches import Rectangle, FancyBboxPatch
    tol_px = 3
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
        leaks = sum(int(x > 0) for x in (ol, orr, ob, ot))
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
    """v3.2 新增：tick 刻度标签（xtick/ytick）两两重叠检测。
    Fig2 panel a 的 x 轴两行式刻度标签（如 "B_naive\n(n=876, OneK1K)"）
    曾因第二行 15 字符超 cell 中心距而互相重叠 34px——check_text 只查
    ax.texts（tick labels 不属其中），故需单独枚举 tick label bbox。
    """
    ticks = list(ax.get_xticklabels()) + list(ax.get_yticklabels())
    ticks = [t for t in ticks if not skip_empty(t)]
    n = len(ticks); over = 0
    for i in range(n):
        ri = ticks[i].get_window_extent()
        for j in range(i + 1, n):
            rj = ticks[j].get_window_extent()
            w, h = inter_wh(ri, rj)
            if w <= 1 or h <= 1:
                continue
            over += 1
            a = ticks[i].get_text().replace("\n", " ")[:30]
            b = ticks[j].get_text().replace("\n", " ")[:30]
            print(f"  TICK {w*h:8.0f}px2 | \"{a}\" <-> \"{b}\"")
    print(f"=== panel {name}: {n} ticklabels, {over} tick-tick overlaps")
    return over


total = 0
for k, ax in enumerate(axs):
    total += check_text(ax, chr(97 + k))
    total += check_patch(ax, chr(97 + k))
    total += check_containment(ax, chr(97 + k))
    total += check_ticklabels(ax, chr(97 + k))
print(f"TOTAL overlaps: {total}")

# ---- panel b 柱颜色像素核对（canvas buffer，避免 bbox=tight 偏移） ----
buf = np.asarray(fig.canvas.buffer_rgba())[:, :, :3]  # (H, W, 3)
leads_col = [
    ("rs13226913", "#0072B2", "C1GALT1"),
    ("rs10238682", "#D55E00", "C1GALT1"),
    ("rs7856182",  "#E69F00", "GALNT12"),
    ("rs5910940",  "#009E73", "C1GALT1C1"),
]
ctxs = ["OneK1K_B_naive", "OneK1K_B_memory", "OneK1K_B_intermediate",
        "GTEx_v10_blood", "GTEx_v10_kidney_cortex"]
bar_w = 0.15; step = 0.24
df = F.pd.read_csv(F.F2B_MATRIX, sep="\t")
axb = axs[1]
mism = 0
print("\n--- panel b bar color pixel check ---")
for j, ctx in enumerate(ctxs):
    for i, (rsid, hexc, gene) in enumerate(leads_col):
        row = df[(df["lead_rsid"] == rsid) & (df["dataset"] == ctx)]
        if row.empty:
            continue
        tested = str(row.iloc[0]["axis_gene_tested"]).strip().lower() == "yes"
        xc = j + (i - 1.5) * step
        if not tested:
            # NA 灰格 = 灰底 #EEEEEE + hatch 线 #999999，采样点可能落在任一
            # 一种上，hash 图样无"真色"，跳过（非视觉问题）
            print(f"  {ctx[:18]:18s} {rsid:12s} NA hatched cell, skip")
            continue
        pv = float(row.iloc[0]["axis_min_p"])
        v = -np.log10(pv) if pv > 0 else 0.0
        if v < 0.5:
            # 柱高 <0.5 unit（约 5-8px）几乎被黑色 edgecolor 覆盖，
            # 像素采样不可靠 → 跳过（v2.0 亦然，非混色）
            print(f"  {ctx[:18]:18s} {rsid:12s} v={v:.2f} too short, skip")
            continue
        yy = v * 0.5
        exp = tuple(int(hexc[k:k+2], 16) for k in (1, 3, 5))
        px, py = axb.transData.transform((xc, yy))
        pxi = int(round(px))
        pyi = buf.shape[0] - 1 - int(round(py))  # buffer 行序 top→bottom，需翻转
        got = tuple(int(x) for x in buf[pyi, pxi]) if (0 <= pxi < buf.shape[1] and 0 <= pyi < buf.shape[0]) else ("OUT",)
        ok = "OK " if (len(got) == 3 and all(abs(g - e) <= 40 for g, e in zip(got, exp))) else "DIFF"
        if ok == "DIFF":
            mism += 1
        print(f"  {ctx[:18]:18s} {rsid:12s} exp={exp} got={got} {ok}")
print(f"panel b color mismatches: {mism}")
print(f"\nALL CHECKED. total overlap={total}, color mismatch={mism}")
