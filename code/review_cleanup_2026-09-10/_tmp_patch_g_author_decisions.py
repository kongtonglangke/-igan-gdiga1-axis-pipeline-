# -*- coding: utf-8 -*-
"""_tmp_patch_g_author_decisions.py — 作者决策 ②③④ 落实（2026-09-10）

② 通讯 Yan Wang ≠ [8] Wang YN（作者确认"否"）→ 无需披露；本条只做记录。
③ GA 图内 "(L7)" 已去除（另见 _tmp_patch_f_ga_delL7.py + 重渲染 GA_v3.3）
   → 本轮同步：ch1 引用 "figure v3.2"→"v3.3"；图注草案 v02→v03（去内部编号+英式拼写）。
④ Additional_File_1_Figure_Legends.md 的 `>` 编辑版本注头块投稿前删
   → 抽出入档 review_cleanup_2026-09-10/ 后从交付件删除。
"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import io, os

ROOT = _paths.LEGACY
PKG = _paths.LEGACY  # 原包根 阶段4.5_复现包
RC = PKG("reproducibility", "code", "review_cleanup_2026-09-10")
GA_DIR = ROOT("阶段4", "主图", "data_check")


def rd(p):
    return io.open(p, encoding="utf-8", newline="").read()


def wr(p, t):
    io.open(p, "w", encoding="utf-8", newline="").write(t)


def rep1(t, old, new, label, require=True):
    n = t.count(old)
    if n != 1:
        print(("  MISS " if require else "  skip ") + f"{label} (count={n})")
        return t, False
    print(f"  OK   {label}")
    return t.replace(old, new, 1), True


# ---------------------------------------------------------------- ① ch1
print("[1] ch1  版本引用 + 修订注")
p1 = PKG("submission", "manuscript", "01_Title_Abstract_Keywords.md")
t1 = rd(p1)

anchor = ("> Revision in v1.6 (2026-09-10, residual internal-code sweep): the word-count note "
          "no longer refers to the internal layer tag (\"the L7 integration\" → "
          "\"the upstream-regulator integration\").")
note = ("\n> Revision in v1.7 (2026-09-10, author decision ③): the Graphical-abstract version "
        "reference was updated from v3.2 to v3.3 — the internal analysis tag \"(L7)\" was removed "
        "from the annotation strip rendered inside the GA artwork, so the submitted graphical "
        "abstract now carries no internal analysis code. No other text changes.")
t1, _ = rep1(t1, anchor, anchor + note, "插入 v1.7 修订注")
t1, _ = rep1(t1, "## Graphical abstract (final design — file GA.pdf; figure v3.2)",
             "## Graphical abstract (final design — file GA.pdf; figure v3.3)",
             "GA 版本 v3.2→v3.3")
wr(p1, t1)

# ---------------------------------------------------------------- ③ GA 图注 v02 → v03
print("\n[2] GA 图注草案 v02 → v03")
p2 = os.path.join(GA_DIR, "ga_caption_v02.md")
t2 = rd(p2)
subs = [
    ("# Graphical Abstract — 图注草案 v02（2026-09-10 · 对应 GA_v3.2）",
     "# Graphical Abstract — 图注草案 v03（2026-09-10 · 对应 GA_v3.3）"),
    ("> **v01 → v02 变更**：GA 图件由 v1.0/v3.1 概念稿升级为 **v3.2（方案 A · L7 上游调控层回补）**——右侧 \"Epigenetic window\" 方框下方新增一枚浅绿注释条\n"
     "> `upstream circuit (L7) / cytokine → DNMT1 → Sp/KLF`，并以绿色箭头指入该方框（语义：回路作用于 C1GALT1C1 启动子，即干预窗所在位点）。\n"
     "> 图注据此补 L7 句；同时把 v01 中\"4 面板 banner / 244×80 mm\"的早期描述更正为实际单幅三块式版面 **244×95 mm**。\n"
     "> v01 留档不删。",
     "> **v02 → v03 变更（作者决策③「去编号」）**：GA 图件由 **v3.2** 升级为 **v3.3**——右侧 \"Epigenetic window\" 方框下方那枚浅绿注释条上的内部层号已去除，\n"
     "> 文字由 `upstream circuit (L7) / cytokine → DNMT1 → Sp/KLF` 改为 `upstream circuit / cytokine → DNMT1 → Sp/KLF`，并以绿色箭头指入该方框（语义：回路作用于 C1GALT1C1 启动子，即干预窗所在位点）。\n"
     "> 图注同步去编号（不再出现 evidence layer L7）；美式 colocalization 统一为英式 colocalisation，与全稿一致。几何/配色/箭头均未变。\n"
     "> **v01 / v02 留档不删。**（v02 的沿革：GA 图件由 v1.0/v3.1 概念稿升级为 v3.2 方案 A · 上游调控层回补——新增该注释条；并把 v01 中\"4 面板 banner / 244×80 mm\"的早期描述更正为实际单幅三块式版面 244×95 mm。）"),
    ("**数据源**：无数据锚点；概念内容 = Background 多击模型 + R3（无共享结构）+ R7（表观干预窗）+ **R8（L7 上游回路）**。",
     "**数据源**：无数据锚点；概念内容 = Background 多击模型 + 无共享因果结构结论（coloc 阴性）+ 获得性表观干预窗 + 上游调控回路（单细胞 regulon 推断）。"),
    ("which sits downstream of an upstream regulatory circuit in which cytokine signalling converges on DNMT1 and Sp/KLF transcription factors at that promoter (single-cell regulon inference; evidence layer L7).",
     "which sits downstream of an upstream regulatory circuit in which cytokine signalling converges on DNMT1 and Sp/KLF transcription factors at that promoter (single-cell regulon inference)."),
    ("(colocalization and allelic-score null)", "(colocalisation and allelic-score null)"),
    ("Lower-right strip: the upstream regulatory circuit (L7) feeding that window.",
     "Lower-right strip: the upstream regulatory circuit feeding that window."),
    ("- **upstream circuit (L7)**：单细胞 CollecTRI regulon 推断所得上游调控层——细胞因子（IL-4/IL-17 等）信号经 DNMT1 与 Sp/KLF 家族转录因子汇聚于 *C1GALT1C1* 启动子（Results R8；Additional file 1: Figures S5–S6、Tables S9–S10）。**该层为相关性证据，不作因果方向声明。**",
     "- **upstream circuit**：单细胞 CollecTRI regulon 推断所得上游调控层——细胞因子（IL-4/IL-17 等）信号经 DNMT1 与 Sp/KLF 家族转录因子汇聚于 *C1GALT1C1* 启动子（Additional file 1: Figures S5–S6、Tables S9–S10）。**该层为相关性证据，不作因果方向声明。**"),
    ("2. **用词**：5-azacytidine / EGCG / no shared causal structure / upstream circuit 与正文 R3/R7/R8 一致；无论文1 数据引用 ✓",
     "2. **用词**：5-azacytidine / EGCG / no shared causal structure / upstream circuit 与正文相应结论小节一致；无论文1 数据引用 ✓"),
    ("3. **L7 表述**：注释条与图注均只写\"上游调控/回路\"的相关性结论，未写\"驱动/因果\"；与正文 R8 及 Limitations 第 6 条的口径一致 ✓",
     "3. **上游回路表述**：注释条与图注均只写\"上游调控/回路\"的相关性结论，未写\"驱动/因果\"；与正文上游调控小节及 Limitations 相关条款的口径一致 ✓"),
    ("> 交付物：`阶段4/主图/out/GA_v3.2.{png,pdf,tiff}`（投稿件已拷为 `submission/figures/GA.{png,pdf,tiff}`）；v3.0/v3.1 留档不删。本图注 v02 即整合用终稿文本。",
     "> 交付物：`阶段4/主图/out/GA_v3.3.{png,pdf,tiff}`（投稿件已拷为 `submission/figures/GA.{png,pdf,tiff}`）；v3.0/v3.1/v3.2 留档不删。本图注 **v03** 即整合用终稿文本。"),
]
for old, new in subs:
    t2, _ = rep1(t2, old, new, old[:46].replace("\n", " "))

# 兜底：确保 v03 内不再有 "L7"
left = t2.count("L7")
print(f"  v03 中残留 'L7' 次数 = {left}（应为 0）")
p3 = os.path.join(GA_DIR, "ga_caption_v03.md")
wr(p3, t2)
print("  wrote:", p3)

# ---------------------------------------------------------------- ④ Additional file 1 头块
print("\n[3] Additional_File_1_Figure_Legends.md 去 `>` 头块（入档后删）")
p4 = PKG("submission", "additional_file_1", "Additional_File_1_Figure_Legends.md")
t4 = rd(p4)
lines = t4.split("\n")
head_idx = [i for i, l in enumerate(lines[:12]) if l.lstrip().startswith(">")]
block = [lines[i] for i in head_idx]
arch = os.path.join(RC, "Additional_File_1_Figure_Legends_header_notes_archived.md")
wr(arch,
   "# [ARCHIVED] Additional file 1 — 图注文件头部编辑版本注（投稿前已从交付件删除）\n\n"
   "> 归档时间：2026-09-10。作者决策④：补充材料**直接交付件**的 `>` 编辑版本注在投稿前删除，\n"
   "> 内容原样留档于此，避免版本信息丢失（母本 `<!-- ASSEMBLY NOTE -->` 同理，见 `_build_submission_fulltext.py`）。\n\n"
   "原文（`submission/additional_file_1/Additional_File_1_Figure_Legends.md` 第 3–7 行）：\n\n"
   + "\n".join(block) + "\n")
print("  archived ->", arch, f"({len(block)} 行)")

keep = [l for l in lines if not (l.lstrip().startswith(">") and lines.index(l) < 12)]
# 去掉因删行产生的连续空行（只处理头部区域）
out, blank = [], 0
for l in keep:
    if l.strip() == "":
        blank += 1
        if blank <= 1:
            out.append(l)
    else:
        blank = 0
        out.append(l)
wr(p4, "\n".join(out))
print("  chk header now:", "\n".join(rd(p4).split("\n")[:6]))
print("\ndone")
