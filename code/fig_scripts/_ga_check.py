"""_ga_check.py — Graphical Abstract 布局自检（text/patch/containment 三件套）

GA 是纯概念示意：单 ax、axis off（无 tick/legend）、3 个 FancyBboxPatch
圆角 box + 6 处 ax.text（其中 3 处 annotate 空文本伪 text）+ 2 条箭头
patch + 1 条分隔虚线。
特殊处理：
  - check_text：跳过 annotate 空文本（get_text().strip() 为空）
  - check_patch：只对 3 个实心 FancyBboxPatch（facecolor alpha=1）两两
    求交，跳过 FancyArrowPatch——箭头端点设计上必须落在 box 边缘/内部
  - check_containment：文字中心落入某实心 box → bbox 必须完全在 box 内
    （GA 文字全部为"文字在盒内"设计，无越盒文字）
  - containment 补充：所有文字 bbox 不得越出画布数据区 [0,10]×[0,4]
    （save_final bbox_inches='tight' 会把越界文字裁进来，说明文字
    跑出设计区 = 布局 bug）
figsize/GridSpec 必须与 ga_render.py main() 一致：
  plt.figure(figsize=(244,95)mm) + add_subplot(111)
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
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import ga_render as G
from style_common import mm_to_in


def iw(r1, r2):
    return min(r1.x1, r2.x1) - max(r1.x0, r2.x0), min(r1.y1, r2.y1) - max(r1.y0, r2.y0)


def check_text(ax, name="GA"):
    over = 0
    all_t = [t for t in ax.texts if t.get_text().strip()]
    for i in range(len(all_t)):
        for j in range(i + 1, len(all_t)):
            r1 = all_t[i].get_window_extent()
            r2 = all_t[j].get_window_extent()
            w, h = iw(r1, r2)
            if w <= 1 or h <= 1:
                continue
            over += 1
            a = all_t[i].get_text().replace("\n", " ")[:30]
            b = all_t[j].get_text().replace("\n", " ")[:30]
            print(f"  TXT {w:5.0f}x{h:3.0f}px | \"{a}\" <-> \"{b}\"")
    print(f"=== {name}: {len(all_t)} texts, {over} text-text overlaps")
    return over


def check_patch(ax, name="GA"):
    """只对不透明实心 box 两两求交；FancyArrowPatch 端点在 box 上属设计意图。"""
    over = 0
    boxes = []
    for p in ax.patches:
        if isinstance(p, FancyArrowPatch):
            continue
        try:
            fc = p.get_facecolor()
        except Exception:
            continue
        if len(fc) == 4 and fc[3] >= 0.95:
            boxes.append(p)
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            r1 = boxes[i].get_window_extent()
            r2 = boxes[j].get_window_extent()
            w, h = iw(r1, r2)
            if w <= 6 or h <= 6:
                continue
            over += 1
            print(f"  PATCH {w:5.0f}x{h:3.0f}px | box#{i} <-> box#{j}")
    print(f"=== {name}: {len(boxes)} solid boxes, {over} box-box overlaps")
    return over


def check_containment(ax, name="GA"):
    """文字中心在实心 box 内 → bbox 必须完全在 box 内。
    另：所有文字 bbox 不得越出画布 [0,10]×[0,4]（bbox_inches='tight' 会把
    越界文字裁回画布，说明文字本不该伸出设计区）。"""
    leaks = 0
    texts = [t for t in ax.texts if t.get_text().strip()]
    boxes = [p for p in ax.patches
             if not isinstance(p, FancyArrowPatch)
             and len(p.get_facecolor()) == 4 and p.get_facecolor()[3] >= 0.95]
    for t in texts:
        r = t.get_window_extent()
        cx, cy = (r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2
        for b in boxes:
            br = b.get_window_extent()
            if br.x0 <= cx <= br.x1 and br.y0 <= cy <= br.y1:
                if (r.x0 < br.x0 - 3 or r.x1 > br.x1 + 3
                        or r.y0 < br.y0 - 3 or r.y1 > br.y1 + 3):
                    leaks += 1
                    print(f"  ESCAPE-BOX {r.width:.0f}x{r.height:.0f}px | "
                          f"\"{t.get_text()[:28].replace(chr(10), ' ')}\" overflow box")
                break  # 中心只可能在一个 box 内
    # 画布越界：把 ax 数据区 [0,10]×[0,4] 换算成窗口坐标
    x0, y0 = ax.transData.transform((0, 0))
    x1, y1 = ax.transData.transform((10, 4))
    x0, x1 = min(x0, x1), max(x0, x1)
    y0, y1 = min(y0, y1), max(y0, y1)
    for t in texts:
        r = t.get_window_extent()
        if (r.x0 < x0 - 3 or r.x1 > x1 + 3
                or r.y0 < y0 - 3 or r.y1 > y1 + 3):
            leaks += 1
            print(f"  ESCAPE-CANVAS {r.width:.0f}x{r.height:.0f}px | "
                  f"\"{t.get_text()[:28].replace(chr(10), ' ')}\" beyond data area")
    print(f"=== {name}: containment leaks={leaks}")
    return leaks


def main():
    fig = plt.figure(figsize=(mm_to_in(244), mm_to_in(95)), dpi=300)
    ax = fig.add_subplot(111)
    G.draw_ga(ax)
    fig.canvas.draw()

    total = 0
    total += check_text(ax)
    total += check_patch(ax)
    total += check_containment(ax)
    print(f"\nTOTAL issues: {total}")


if __name__ == "__main__":
    main()
