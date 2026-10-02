# -*- coding: utf-8 -*-
"""
_tmp_patch_t_round3b_author_decisions.py
第三轮核查「作者裁定」办结（2026-09-10 深夜 · 包 v2.4）

作者两项裁定：
  ① 单细胞病例对照实算（报告 §1.2）→ **不写**入稿件。仅留档于报告与本文 Reminder，稿件零改动。
  ② Discussion 末段「实验室在研正交实验」前瞻段 → **删除**。
     裁定理由：作者实验室**并未**开展该组实验，原句属不实陈述。
     同时按作者意见把该段收尾改写为不突兀的段落结语。

改动范围（唯一）：submission/manuscript/05_Discussion_Limitations.md
  · 段 5（Discussion 末段）尾句由「we state them as explicit tests for the orthogonal
    systems under way in our laboratory: (i)…(ii)…(iii)… and will be reported separately.」
    改写为「不主张本组在研、而把三个预测作为该框架可被检验的方向 + 收束全段的结语」。
  · 章节头补 v1.8 版本注。
  数值、判定、图件、图注、结论零改动。

用法：python reproducibility/code/review_cleanup_2026-09-10/_tmp_patch_t_round3b_author_decisions.py
"""
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))  # 阶段4.5_复现包
M = os.path.join(BASE, "submission", "manuscript")
CH5 = os.path.join(M, "05_Discussion_Limitations.md")

START = "The logical next steps are therefore functional, and we state them as"
NEW = (
    "The logical next steps are therefore functional. The circuit defined here makes specific, "
    "falsifiable predictions\u2014that restoring C1GALT1C1 expression, or reversing promoter "
    "hypermethylation, restores galactosylation capacity in B cells from patients with IgAN, and that "
    "the predicted methylation-sensitive KLF2 and SP1 motifs at the C1GALT1C1 promoter are "
    "occupancy-blocked by methylation\u2014that perturbation and chromatin-level experiments can now "
    "test directly. Closing that gap is the one thing a design built entirely on public summary data "
    "cannot deliver, and it defines the natural direction for the next studies of the axis. What this "
    "framework contributes is therefore a boundary and a direction: it puts a boundary on how much of "
    "IgAN susceptibility common germline variation in the axis can explain, and it identifies the "
    "acquired epigenetic state of the B-cell compartment as the interface at which the axis becomes "
    "both experimentally testable and therapeutically addressable."
)

VN_MARK = "> v1.7 (2026-09-10, third reviewer-perspective pass)"
VN_NEW = (
    "> v1.8 (2026-09-10, author decision on the closing paragraph): the final paragraph no longer "
    "states that the three functional tests are \"under way in our laboratory\" or that they \"will be "
    "reported separately\"; the experimental predictions are now presented as what the framework makes "
    "directly testable and as the direction for the next studies of the axis, and the paragraph closes "
    "on the boundary-and-direction summary of the study's contribution. The change removes a claim about "
    "the authors' own in-progress work that could not be substantiated; no data, values or conclusions changed."
)


def main():
    s = open(CH5, encoding="utf-8").read()

    # --- 1. 尾段替换（定位到行尾，避免手抄 em dash 出错；幂等）---
    hits = [m.start() for m in re.finditer(re.escape(START), s)]
    if len(hits) == 1:
        i = hits[0]
        j = s.find("\n", i)
        if j == -1:
            j = len(s)
        old = s[i:j]
        for must in ("under way in our laboratory", "will be reported separately"):
            if must not in old:
                print("FAIL: expected phrase missing from the slice:", must); sys.exit(1)
        s = s[:i] + NEW + s[j:]
        print("step1: tail paragraph replaced")
    elif NEW in s:
        print("step1: already applied (idempotent skip)")
    else:
        print("FAIL: neither the original marker nor the new text found"); sys.exit(1)

    # --- 2. 章节头 v1.8 版本注（幂等：已存在则跳过）---
    lines = s.split("\n")
    if not any(l.startswith("> v1.8 (2026-09-10, author decision") for l in lines):
        idx = [k for k, l in enumerate(lines) if l.startswith(VN_MARK)]
        if len(idx) != 1:
            print("FAIL: v1.7 version-note line hit count =", len(idx)); sys.exit(1)
        lines.insert(idx[0] + 1, VN_NEW)
        s = "\n".join(lines)
        print("step2: v1.8 version note inserted")
    else:
        print("step2: v1.8 version note already present (idempotent skip)")

    open(CH5, "w", encoding="utf-8", newline="\n").write(s)
    print("patched:", CH5)

    # --- 3. 复核：旧串在 ch5 正文（剥离 > 版本注行）内清零 ---
    body_only = "\n".join(l for l in open(CH5, encoding="utf-8").read().split("\n")
                          if not l.startswith(">"))
    for k in ("under way in our laboratory", "will be reported separately",
              "explicit tests for the orthogonal"):
        n = body_only.count(k)
        print(f"  stale in ch5 body {k!r}: {n} -> {'PASS' if n == 0 else 'FAIL'}")
        if n:
            sys.exit(1)

    # --- 4. 重建母本 ---
    r = subprocess.run([sys.executable, os.path.join(BASE, "reproducibility", "code",
                                                     "_build_submission_fulltext.py")],
                       cwd=BASE, capture_output=True, text=True)
    print(r.stdout.strip() or r.stderr.strip())

    # --- 5. 词数复算（口径＝原始空白分词；排除头部 > 注 / --- 分隔线）---
    def wc(seg, drop_heading=False):
        ls = [l for l in seg.split("\n") if l.strip()]
        if drop_heading:
            ls = [l for l in ls if not l.strip().startswith("#")]
        return len(" ".join(ls).split())

    body = [l for l in open(CH5, encoding="utf-8").read().split("\n")
            if not l.startswith(">") and not re.fullmatch(r"\s*-{3,}\s*", l)]
    t = "\n".join(body)
    a = t.find("## Discussion"); b = t.find("## Limitations")
    disc, lim = t[a:b], t[b:]
    print(f"ch5 word counts: Discussion {wc(disc)} incl / {wc(disc, True)} excl heading; "
          f"Limitations {wc(lim)} / {wc(lim, True)}")


if __name__ == "__main__":
    main()
