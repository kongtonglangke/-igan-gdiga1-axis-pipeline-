#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""_tmp_patch_e_residual_Rn.py — 补漏：清除上一轮遗漏的「裸 Rn」内部编号
上一轮只清了括号式 (R2)，漏了 'reported in R2' / 'Corresponds to R2' 这类裸写法。
本轮仅改投稿正文内真泄漏（ch3/ch4 正文各 1 处、ch8 表注 6 处）+ ch1 一处非投稿字数注。
内部集成附录（ch4 的 'Figure calls in text' 行与 reproducibility/data/M*_ 路径）保留不动。
"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import os

MS = _paths.at("阶段4.5_复现包/submission/manuscript")


def patch(fn, pairs):
    p = MS(fn)
    with open(p, encoding="utf-8", newline="") as fh:
        t = fh.read()
    for old, new, exp in pairs:
        n = t.count(old)
        if n != exp:
            raise SystemExit(f"[FAIL] {fn}: found {n}, expected {exp}\n  >> {old[:130]}")
        t = t.replace(old, new)
    with open(p, "w", encoding="utf-8", newline="") as fh:
        fh.write(t)
    print(f"OK  {fn}  ({len(pairs)} pairs)")


# ============================== ch1 ==============================
patch("01_Title_Abstract_Keywords.md", [
    ("were re-measured after the L7 integration and after the 2026-09-10 reviewer-perspective clean-up",
     "were re-measured after the upstream-regulator integration and after the 2026-09-10 reviewer-perspective clean-up", 1),
    ("> Revision in v1.5 (2026-09-10, reviewer-perspective clean-up):",
     "> Revision in v1.6 (2026-09-10, residual internal-code sweep): the word-count note no longer refers to the internal layer tag (\"the L7 integration\" \u2192 \"the upstream-regulator integration\").\n> Revision in v1.5 (2026-09-10, reviewer-perspective clean-up):", 1),
])

# ============================== ch3 ==============================
patch("03_Methods.md", [
    ("independent corroboration of the cell-state-specific eQTL signal reported in R2; it is not a quantitative estimate",
     "independent corroboration of the cell-state-specific eQTL signal reported above; it is not a quantitative estimate", 1),
    ("> v1.8 (2026-09-10, reviewer-perspective clean-up):",
     "> v1.9 (2026-09-10, residual internal-code sweep): the closing sentence of the single-cell subsection no longer carries a bare result tag (\"reported in R2\" \u2192 \"reported above\").\n> v1.8 (2026-09-10, reviewer-perspective clean-up):", 1),
])

# ============================== ch4 ==============================
patch("04_Results.md", [
    ("(the marker-gated partition of R2 assigns more cells to the plasma-like compartment",
     "(the marker-gated partition assigns more cells to the plasma-like compartment", 1),
    ("> v1.7 (2026-09-10, reviewer-perspective clean-up):",
     "> v1.8 (2026-09-10, residual internal-code sweep): the single-cell atlas paragraph no longer carries a bare result tag (\"the marker-gated partition of R2\" \u2192 \"the marker-gated partition\"). The internal \"Figure calls in text\" appendix line and the reproducibility-package paths (reproducibility/data/M1_geo ... M6_intervention) intentionally retain the original analysis-directory names, because the code and the reproducibility README resolve those paths; they are stripped from the assembled submission full text.\n> v1.7 (2026-09-10, reviewer-perspective clean-up):", 1),
])

# ============================== ch8 ==============================
patch("08_Figure_Legends_Tables.md", [
    ("Corresponds to R1 and Fig. 1b.", "Corresponds to Fig. 1b.", 1),
    ("Corresponds to R2 and Fig. 2b.", "Corresponds to Fig. 2b.", 1),
    ("Corresponds to R3 and Fig. 3b, c.", "Corresponds to Fig. 3b, c.", 1),
    ("Corresponds to R4 and Fig. 4a\u2013c.", "Corresponds to Fig. 4a\u2013c.", 1),
    ("Corresponds to R6 and Fig. 5a.", "Corresponds to Fig. 5a.", 1),
    ("Corresponds to R7 and Fig. 6a.", "Corresponds to Fig. 6a.", 1),
    ("> v1.7 (2026-09-10, L7 null re-run at 1,000 sets):",
     "> v1.8 (2026-09-10, residual internal-code sweep): the six main-text table notes now read \"Corresponds to Fig. X.\" without the internal result tag (\"Corresponds to R1 and Fig. 1b\" \u2192 \"Corresponds to Fig. 1b\", and likewise for R2/R3/R4/R6/R7); the Fig. 1c legend keeps the in-panel layer index (L0\u2013L7), which the panel itself prints.\n> v1.7 (2026-09-10, L7 null re-run at 1,000 sets):", 1),
])

print("done")
