#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""_readme_update_v1.3.py — reproducibility/README 升 v1.3"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import os

ROOT = _paths.at("阶段4.5_复现包/reproducibility")


def patch(fn, pairs):
    p = ROOT(fn)
    with open(p, encoding="utf-8", newline="") as fh:
        t = fh.read()
    for old, new, exp in pairs:
        n = t.count(old)
        if n != exp:
            raise SystemExit(f"[FAIL] {fn}: found {n}, expected {exp}\n  >> {old[:150]}")
        t = t.replace(old, new)
    with open(p, "w", encoding="utf-8", newline="") as fh:
        fh.write(t)
    print(f"OK  {fn}  ({len(pairs)} pairs)")


patch("README.md", [
    ("**审稿整改补丁脚本 5 件**（体例统一、去内部编号、定点改写、L7 数字同步、Additional file 图注修正）",
     "**审稿整改补丁脚本 6 件**（体例统一、去内部编号、定点改写、L7 数字同步、Additional file 图注修正、**裸 Rn 补漏**）", 1),

    ("说明内部编号仅存于本文件的定位口径。",
     "说明内部编号仅存于本文件的定位口径。",
     1),

    ("*v1.2 · 2026-09-10 · 阶段4.5 复现包* — v1.1→v1.2：补列 `sc_regulatory/`、`directional_validation/`、`C1_eas_eqtl/`、`review_cleanup_2026-09-10/` 四个子目录；登记 step2 `--k-null 1000` 重跑与 motif 交集因子变化；说明内部编号仅存于本文件的定位口径。",
     "*v1.3 · 2026-09-10 · 阶段4.5 复现包* — v1.2→v1.3：补登 `review_cleanup_2026-09-10/` 第 6 件补丁脚本（`_tmp_patch_e_residual_Rn.py`，清除上一轮遗漏的 8 处裸 `Rn`）。另说明两处**有意保留**的内部编号（**均不进母本**、拼装时自动剔除）：ch4 集成附录的 \u201cFigure calls in text\u201d 行，以及 `data/M1_geo … M6_intervention` 数据目录路径（代码与本 README 依赖其原名，重命名会断脚本）。\n"
     "*v1.2 · 2026-09-10 · 阶段4.5 复现包* — v1.1→v1.2：补列 `sc_regulatory/`、`directional_validation/`、`C1_eas_eqtl/`、`review_cleanup_2026-09-10/` 四个子目录；登记 step2 `--k-null 1000` 重跑与 motif 交集因子变化；说明内部编号仅存于本文件的定位口径。", 1),
])

print("done reproducibility/README.md")
