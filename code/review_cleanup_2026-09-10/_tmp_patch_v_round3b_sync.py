# -*- coding: utf-8 -*-
"""
_tmp_patch_v_round3b_sync.py
第三轮作者裁定办结（② 删「实验室在研实验」前瞻段）——包内文档与报告同步 + 临时件清理。

改动：
  · README.md（包）            v2.3 → v2.4：标题、组件快照（ch5 v1.8 / Discussion 1,600·1,598 /
                               正文合计 9,952·9,722）、目录说明计数、全稿版本列表、遗留清理、版本记录新增 v2.4 条。
  · reproducibility/README.md  v1.8 → v1.9：装配说明版本、目录说明计数、版本记录新增 v1.9 条。
  · 审稿视角第三轮核查报告_v1.0_2026-09-10.md：追加「§六 作者裁定办结」。
  · 删除包根临时快照 `_ft_before.md`（可由装配脚本随时重建）。

用法：python reproducibility/code/review_cleanup_2026-09-10/_tmp_patch_v_round3b_sync.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))  # 阶段4.5_复现包
PKG_README = os.path.join(BASE, "README.md")
REP_README = os.path.join(BASE, "reproducibility", "README.md")
REPORT = os.path.join(BASE, "审稿视角第三轮核查报告_v1.0_2026-09-10.md")
TEMP_FT = os.path.join(BASE, "_ft_before.md")


def sub(path, pairs, label):
    s = open(path, encoding="utf-8").read()
    for old, new in pairs:
        n = s.count(old)
        if n != 1:
            print(f"FAIL[{label}] hit={n} for: {old[:70]!r}")
            sys.exit(1)
        s = s.replace(old, new)
    open(path, "w", encoding="utf-8", newline="\n").write(s)
    print(f"OK [{label}] {len(pairs)} edit(s) -> {os.path.basename(path)}")


# ---------------- 1. 包 README ----------------
PKG = [
    ("# 阶段 4.5 复现包 / 投稿材料包（v2.3 · 2026-09-10）",
     "# 阶段 4.5 复现包 / 投稿材料包（v2.4 · 2026-09-10）"),

    ("审稿整改补丁脚本 9 件 + README 升版脚本 4 件 + 图注头块归档件 1 件（可审计，非分析流程）",
     "整改留痕脚本（`_tmp_*`）24 件 + 渲染/升版辅助脚本 5 件 + 图注头块归档件 1 件（可审计，非分析流程）"),

    ("| Discussion + Limitations | ch5 **v1.7**（第三轮：去「单调梯度」表述 + 单细胞表达量口径；此前去编号 + causal-inference 边界声明 + **三条验证清单**） | Discussion **1,589/1,587** / Limitations **605/603** | 正文 |",
     "| Discussion + Limitations | ch5 **v1.8**（作者裁定：删去「实验室在研正交实验」前瞻段（该组实验实未开展）、末段收尾改写为「边界与方向」结语；第三轮：去「单调梯度」表述 + 单细胞表达量口径；此前去编号 + causal-inference 边界声明 + 三条验证清单） | Discussion **1,600/1,598** / Limitations **605/603** | 正文 |"),

    ("| **正文合计（章节2–6，含小节标题）** | — | **9,941 含 / 9,711 不含**（= BG 762/760 + Methods 3,626/3,521 + Results 3,211/3,094 + Disc 1,589/1,587 + Lim 605/603 + Concl 148/146；2026-09-10 第三轮复核后按各章增量实测） | 脚本实测定格 |",
     "| **正文合计（章节2–6，含小节标题）** | — | **9,952 含 / 9,722 不含**（= BG 762/760 + Methods 3,626/3,521 + Results 3,211/3,094 + Disc 1,600/1,598 + Lim 605/603 + Concl 148/146；2026-09-10 第三轮复核 + 作者裁定收尾改写后按各章增量实测） | 脚本实测定格 |"),

    ("附：本稿正文合计 **9,941** 含小节标题；",
     "附：本稿正文合计 **9,952** 含小节标题；"),

    ("v1.6（2026-09-10）：**补漏重建**——清除上一轮遗漏的 8 处**裸 `Rn`**（ch3/ch4 正文各 1 处、ch8 表注 6 处）；重建后正文 109,804 chars，母本裸 `Rn`/`Mn` = 0。",
     "v1.6（2026-09-10）：**补漏重建**——清除上一轮遗漏的 8 处**裸 `Rn`**（ch3/ch4 正文各 1 处、ch8 表注 6 处）；重建后正文 109,804 chars，母本裸 `Rn`/`Mn` = 0。\n"
     "v1.7 / v1.8 / v1.9 / v1.10（2026-09-10）：装配说明随各轮整改同步升版——图表发表性审查（v1.7，母本 **110,787** chars）、作者三项裁定 ①改图注 ②八层 ③补充表清理（v1.8，**111,565** chars）、第三轮审稿视角十处正文整改（v1.9，**112,900** chars）、作者裁定「删实验室在研实验前瞻段」（v1.10，**113,499** chars，其中新增装配说明注约 555 chars、不进稿件正文；本轮正文净增 **+44** chars）。"),

    ("- **遗留清理（可选）**：包根 `_ft_before.md`（旧母本快照）与 `reproducibility/code/_tmp_build_tables.py`（ch8 Tables 转录脚本，README 已登记）投稿前可择一处理——前者建议删除、后者建议保留并去掉 `_tmp_` 前缀以明示其为正式转录脚本。",
     "- **遗留清理**：包根 `_ft_before.md`（旧母本快照）**已于 2026-09-10 深夜删除**（属临时快照；母本随时可由 `reproducibility/code/_build_submission_fulltext.py` 重建）；`reproducibility/code/_tmp_build_tables.py`（ch8 Tables 转录脚本，README 已登记）建议保留并去掉 `_tmp_` 前缀以明示其为正式转录脚本。"),
]

V24 = (
    "- v2.4（2026-09-10 深夜）：**第三轮核查「作者裁定」办结**。① 裁定一（单细胞病例对照实算是否写入稿件）＝**不写**：报告 §1.2 的 GSE285335 分组检验（阴性、且该库病期成块编号致混杂）仅留档，稿件零改动。"
    "② 裁定二＝**删除 Discussion 末段「实验室在研正交实验」前瞻段**（作者确认该组实验并**未**开展，原「explicit tests … under way in our laboratory」「will be reported separately」属不实陈述）："
    "ch5 末段尾句改写为——三个预测作为该框架**可直接检验**的方向呈现（不再主张本组在研），并以「What this framework contributes is therefore a boundary and a direction: …」收束全段（作者要求结尾不生硬突兀）。"
    "版本：ch5 **v1.7 → v1.8**；Discussion **1,589/1,587 → 1,600/1,598**（+11 含 / +11 不含）；正文合计 **9,941/9,711 → 9,952/9,722**；母本重建 **113,499** chars（装配说明 **v1.9 → v1.10**）。"
    "机检全过：Abstract **348**（限 350）、主文表 1–7 全部被引、引用首现序 **[1..27] 严格递增**、正文中文 **0**、内部编号残留 **0**（仅存 Fig. 1c 图内索引）、陈旧串（含本轮三条：under way in our laboratory / will be reported separately / explicit tests for the orthogonal）**全 0**。"
    "数值、图件、结论零改动。留痕脚本 `_tmp_patch_t_round3b_author_decisions.py`、`_tmp_patch_u_round3b_header_verify.py`、`_tmp_patch_v_round3b_sync.py`；报告 `审稿视角第三轮核查报告_v1.0_2026-09-10.md` 追加 §六。"
)


def patch_pkg():
    sub(PKG_README, PKG, "pkg-readme")
    # 版本记录追加 v2.4（插在 v2.3 条目之后）
    s = open(PKG_README, encoding="utf-8").read()
    anchor = "报告 `审稿视角第三轮核查报告_v1.0_2026-09-10.md`。\n"
    if s.count(anchor) != 1:
        print("FAIL: v2.3 entry anchor hit =", s.count(anchor)); sys.exit(1)
    s = s.replace(anchor, anchor + V24)
    open(PKG_README, "w", encoding="utf-8", newline="\n").write(s)
    print("OK [pkg-readme] v2.4 version-record entry appended")


# ---------------- 2. reproducibility/README ----------------
REP = [
    ("| `code/_build_submission_fulltext.py` | 章节定稿 → 单文件全稿（Manuscript）拼装脚本（头部装配说明 v1.9） | 全稿组装 |",
     "| `code/_build_submission_fulltext.py` | 章节定稿 → 单文件全稿（Manuscript）拼装脚本（头部装配说明 v1.10） | 全稿组装 |"),

    ("| `code/review_cleanup_2026-09-10/` | **审稿整改补丁脚本 9 件**（体例统一、去内部编号、定点改写、L7 数字同步、Additional file 图注修正、**裸 Rn 补漏**、**GA 去编号**、**作者决策 ②③④ 落实**、**补充材料 -isation**）+ README 升版脚本 4 件 + `Additional_File_1_Figure_Legends_header_notes_archiv",
     "| `code/review_cleanup_2026-09-10/` | **整改留痕脚本（`_tmp_*`）24 件**（体例统一、去内部编号、定点改写、L7 数字同步、Additional file 图注修正、裸 Rn 补漏、GA 去编号、作者决策 ②③④ 落实、补充材料 -isation、**第三轮十处整改**、**本轮删实验室在研前瞻段**）+ 渲染/升版辅助脚本 5 件 + `Additional_File_1_Figure_Legends_header_notes_archiv"),
]

V19 = (
    "*v1.9 · 2026-09-10（深夜）· 阶段4.5 复现包* — v1.8→v1.9：**第三轮作者裁定办结的包内同步**。"
    "`_build_submission_fulltext.py` 装配说明 **v1.9→v1.10**（记本轮改动：Discussion 末段删去「under way in our laboratory」「will be reported separately」并改写收尾为「boundary and direction」结语；单细胞病例对照实算经作者裁定不写入稿件），"
    "母本重建 **113,499 chars**；新增留痕脚本 `review_cleanup_2026-09-10/_tmp_patch_t_round3b_author_decisions.py`（ch5 尾段定点改写 + v1.8 版本注，幂等、逐条断言）、`_tmp_patch_u_round3b_header_verify.py`（装配说明升版 + 重建 + 机检 A–H）、`_tmp_patch_v_round3b_sync.py`（本文件与包 README 同步 + 临时件清理）。"
    "删除包根临时快照 `_ft_before.md`（可由装配脚本重建）。**数值、图件、结论零改动。**"
)


def patch_rep():
    sub(REP_README, REP, "rep-readme")
    s = open(REP_README, encoding="utf-8").read()
    anchor = "*v1.8 · 2026-09-10（深夜）· 阶段4.5 复现包*"
    if s.count(anchor) != 1:
        print("FAIL: v1.8 entry anchor hit =", s.count(anchor)); sys.exit(1)
    s = s.replace(anchor, V19 + "\n" + anchor)
    open(REP_README, "w", encoding="utf-8", newline="\n").write(s)
    print("OK [rep-readme] v1.9 version-record entry inserted")


# ---------------- 3. 报告追加 §六 ----------------
SEC6 = """
---

## 六、作者裁定办结（2026-09-10 深夜 · 包 v2.4）

第五节的余留项 ① 与 ② 已由作者裁定并办结：

| # | 裁定 | 处置 |
|---|---|---|
| ① | 单细胞病例对照实算（§1.2）**不写**入稿件 | 稿件**零改动**；该组数字（B 细胞五轴基因 IgAN vs HC 中位 log2 CPM 差 −0.16～+0.46、P = 0.16～0.91，全阴）仅留档于本报告 §1.2 与包内日志，备审稿问答引用 |
| ② | Discussion 末段「实验室在研正交实验」前瞻段 **删除** | 作者确认该组实验**并未开展**，原句「explicit tests … under way in our laboratory」与「will be reported separately」属不实陈述 → ch5 末段尾句改写（见下）|

**②的改写要点（ch5 v1.7 → v1.8）**：删除「本组在研」的归属与「将另行报道」的承诺，把三个实验预测（去甲基化/DNMT1 抑制剂恢复 C1GALT1C1 表达与半乳糖化能力；IL-4/IL-17 驱动的启动子高甲基化可逆性；KLF2/SP1 motif 的甲基化占位阻断）改为**该框架可直接检验的方向**；并按作者「结尾不要生硬突兀」的要求，用一句收束全段：「*What this framework contributes is therefore a boundary and a direction: it puts a boundary on how much of IgAN susceptibility common germline variation in the axis can explain, and it identifies the acquired epigenetic state of the B-cell compartment as the interface at which the axis becomes both experimentally testable and therapeutically addressable.*」

**版本与词数**：ch5 **v1.7 → v1.8**；Discussion **1,589/1,587 → 1,600/1,598**（+11）；Limitations **605/603 不变**；正文合计 **9,941/9,711 → 9,952/9,722**；母本重建 **113,499** chars（装配说明 **v1.9 → v1.10**）。

**机检（`_tmp_patch_u_round3b_header_verify.py`）**：Abstract **348**（≤350）、主文表 1–7 全部被引、引用首现序 **[1..27] 严格递增**、正文中文 **0**、内部编号残留 **0**（仅 Fig. 1c 图内索引）、陈旧串全 0（新增本轮三条：`under way in our laboratory` / `will be reported separately` / `explicit tests for the orthogonal`）→ **ALL CHECKS PASS**。

**余留（更新）**：③ Funding（2027 省基金）④ [26] Fung 投稿当日复核 ⑤ 全体作者终校（署名/CRediT/Competing/Acknowledgements）⑥ GitHub 建仓 → Zenodo DOI 回填 ⑦ .docx/.pdf 生成。

"""


def patch_report():
    s = open(REPORT, encoding="utf-8").read()
    anchor = "*第三轮核查报告 v1.0 · 2026-09-10 深夜 · 数值、图件、结论零改动；改动限于表述准确性、引用完整性与投稿件口径一致性。*"
    if s.count(anchor) != 1:
        print("FAIL: report closing anchor hit =", s.count(anchor)); sys.exit(1)
    s = s.replace(anchor, SEC6.strip("\n") + "\n\n" + anchor)
    # 第五节余留表首列标注办结
    old1 = "| ① | **§1.2 的单细胞病例对照实算是否写入稿件** | 现为不新增"
    new1 = "| ① | ~~**§1.2 的单细胞病例对照实算是否写入稿件**~~ **已办结（不写）** | 现为不新增"
    old2 = "| ② | **Discussion 末段\"explicit tests … under way in our laboratory\"** | 属对作者实验室在研工作的前瞻陈述，无法由公开材料核验，**请作者确认表述与实际相符** |"
    new2 = "| ② | ~~**Discussion 末段\"explicit tests … under way in our laboratory\"**~~ **已办结（删除）** | 作者确认该组实验未开展 → 该句及「will be reported separately」**已删除**，末段改写并加结语（详 §六）|"
    for o, n in ((old1, new1), (old2, new2)):
        if s.count(o) != 1:
            print("FAIL: §5 row anchor hit =", s.count(o), repr(o[:50])); sys.exit(1)
        s = s.replace(o, n)
    open(REPORT, "w", encoding="utf-8", newline="\n").write(s)
    print("OK [report] §六 appended + §5 rows marked closed")


# ---------------- 4. 临时件清理 ----------------
def drop_temp():
    if os.path.exists(TEMP_FT):
        os.remove(TEMP_FT)
        print("OK [cleanup] removed", TEMP_FT)
    else:
        print("OK [cleanup] _ft_before.md already absent")


if __name__ == "__main__":
    patch_pkg()
    patch_rep()
    patch_report()
    drop_temp()
