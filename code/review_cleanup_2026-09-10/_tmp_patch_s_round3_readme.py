# -*- coding: utf-8 -*-
"""
_tmp_patch_s_round3_readme.py — 第三轮整改的包内文档同步

① 阶段4.5_复现包/README.md：标题版本 v2.2→v2.3；刷新「稿件定稿版本快照」表中已落后
   一轮的版本号与词数（ch1/ch2/ch3/ch4/ch5/ch6/ch8、正文合计）；追加 v2.3 版本历史条目。
② reproducibility/README.md：新增 v1.8 条目；`_build_submission_fulltext.py` 的装配说明
   版本 v1.7→v1.9。
每条替换均断言唯一命中。
"""
import os
import sys

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
R1 = os.path.join(BASE, "README.md")
R2 = os.path.join(BASE, "reproducibility", "README.md")

V23 = ("- v2.3（2026-09-10 深夜）：**审稿视角第三轮核查 + 全文整改**。先按「审稿人视角是否还能抬高文章层次」评估："
       "把唯一可行的候选新分析——用已入稿的单细胞数据集（GSE285335，26 供体）做**独立的分组（IgAN vs 健康对照）轴基因表达检验**，以补强"
       "「PBMC 仅 2 例对照」这一最弱环节——实算了一遍（B 细胞 / 浆细胞逐供体 pseudobulk，log2 CPM 与检出率，Mann–Whitney）："
       "结果为**阴性且功效不足**（B 细胞内五个轴基因 IgAN vs HC 中位 log2 CPM 差 −0.16～+0.46，P = 0.16～0.91；唯一 P < 0.05 的是 NK 的 GALNT2，"
       "属非 B 谱系），且该数据集的 26 个样本按疾病分期成块连续编号（GSM8700986–991 晚期 / 992–1002 早期 / 1003–1011 健康），病例对照与处理顺序部分混杂。"
       "该结果**与稿件 Limitations 第七点已自陈的「该数据集病例对照在若干轴基因上与 bulk PBMC 结果方向不一致、故仅用于细胞状态与上游调控刻画」完全一致**，"
       "故**不新增该分析**（新增即是负面而非抬升）；其余可想到的抬升路径（Gd-IgA1 全效应量 MR / coloc、独立 IgAN GWAS 复现、东亚 B 细胞 eQTL 精细定位）"
       "均受「零受控全公开」数据约束所限，不可行。据此转入全维度复核并整改 10 处：**① Table 5 此前未被任何正文句引用**（BMC/GM 要求正文表按序在文内被引用）"
       "→ 甲基化节补 `(Fig. 4a; Table 5)`；**② 同一事实数量级口径矛盾**（ch3「three orders of magnitude」vs ch4「two orders」，实为 9.4e-6 vs 5e-8 ≈ 2.3 个数量级）"
       "→ 两处统一为「approximately two orders of magnitude above the genome-wide significance threshold (P = 9.4 × 10⁻⁶ versus 5 × 10⁻⁸)」；"
       "**③ Discussion 与已修正的 Results/Fig. 2c 冲突**（仍写「attenuate along the naive-to-memory transition」，而图注/正文明确排除单调梯度）"
       "→ 改为「attenuated in the antigen-experienced states, with partial recovery in memory」；**④ Abstract 悬垂分词**（「Produced by …, whether … remains unresolved」）"
       "→ 改为同位语结构，**词数净变 0（348，限 350）**；**⑤ 单细胞表达量口径自相矛盾**（「essentially flat」与「P = 2.9 × 10⁻⁷⁵」并置）→ 改为给出实测均值（1.18 → 1.05、1.02 → 0.85）并说明大细胞数下小位移亦显著；"
       "**⑥ 异质性检验补明所测 lead**（rs13226913，Table S11 中两 lead 各有其值）；**⑦ 缩写表 PPH3/PPH4 → PPH0–PPH4**（与 Results 定义记法一致）；"
       "**⑧ Sp/KLF、IDG、Pharos 四级（Tclin/Tchem/Tbio/Tdark）、NBDC 首次出现处补全称**；**⑨ 缩写表注**关于「所有基因符号均于首用处展开」的表述与实际不符 → 改为准确表述；"
       "**⑩ cover letter 仍写「seven sequential evidence layers」且枚举 (i)–(vii) 漏第 8 层**（正文已统一八层）→ 升 v1.2，改 eight 并补入「疾病层面直接查证」，"
       "疾病阴性证据清单口径写准（colocalisation + 直接回查 + 等位评分 + CpG 甲基化 MR）。版本：ch1 **v1.9** / ch3 **v1.10** / ch4 **v1.9** / ch5 **v1.7** / ch6 **v1.6**"
       "（ch2/ch9 未动）；母本重建 **112,900 chars**（装配说明 **v1.9**）。机检全过：Abstract **348**（限 350）、主文表 1–7 全部被引且按序、引用首现序 **[1..27] 严格递增**、"
       "正文中文 **0**、内部编号残留 **0**（仅存 Fig. 1c 图内索引）、陈旧串全 0。数值、图件、结论零改动。留痕脚本 `_tmp_patch_p_round3_fixes.py`、"
       "`_tmp_patch_q_round3_versions.py`、`_tmp_patch_r_round3_verify.py`；报告 `审稿视角第三轮核查报告_v1.0_2026-09-10.md`。\n")

R2_18 = ("*v1.8 · 2026-09-10（深夜）· 阶段4.5 复现包* — v1.7→v1.8：**审稿视角第三轮整改的包内同步**。"
         "`_build_submission_fulltext.py` 装配说明 **v1.8→v1.9**（记十处正文改动：Table 5 补正文引用、数量级口径统一、Discussion 去「单调梯度」表述、Abstract 悬垂分词修正、"
         "单细胞表达量口径、异质性检验补明 lead、PPH0–PPH4 记法、Sp/KLF 与 IDG/Pharos 四级/NBDC 全称、缩写表注准确性、cover letter 七层→八层），母本重建 **112,900 chars**；"
         "新增留痕脚本 `review_cleanup_2026-09-10/_tmp_patch_p_round3_fixes.py`（18 处定点改写，逐条断言唯一命中）、`_tmp_patch_q_round3_versions.py`（版本注落位）、"
         "`_tmp_patch_r_round3_verify.py`（机检：Abstract 词数 / 主文表正文引用 / 引用首现序 / 中文 / 内部编号 / 陈旧串）。**数值、图件、结论零改动。**\n")

EDITS = [
    (R1, "（v2.2 · 2026-09-10）", "（v2.3 · 2026-09-10）"),
    (R1, "| Title/Abstract/Keywords | ch1 **v1.7**（审稿清洗", "| Title/Abstract/Keywords | ch1 **v1.9**（第三轮复核：Abstract 悬垂分词修正，词数不变；此前审稿清洗"),
    (R1, "| Abstract **346**（限 350 ✅） |", "| Abstract **348**（限 350 ✅） |"),
    (R1, "| Background | ch2 **v1.5**（", "| Background | ch2 **v1.6**（"),
    (R1, "| Methods（含 Table 1） | ch3 **v1.9**（16 小节；", "| Methods（含 Table 1） | ch3 **v1.10**（16 小节；第三轮：数量级口径统一 + IDG/Pharos 四级/NBDC 补全称；"),
    (R1, "| **3,592 含 / 3,487 不含** | 正文 |", "| **3,626 含 / 3,521 不含** | 正文 |"),
    (R1, "| Results R1–R8 | ch4 **v1.8**（去编号", "| Results R1–R8 | ch4 **v1.9**（第三轮：补 Table 5 正文引用 + 异质性检验补明 lead + Sp/KLF 全称；此前去编号"),
    (R1, "| **3,201 含 / 3,084 不含**（8 小节） | 正文 |", "| **3,211 含 / 3,094 不含**（8 小节） | 正文 |"),
    (R1, "| Discussion + Limitations | ch5 **v1.6**（去编号", "| Discussion + Limitations | ch5 **v1.7**（第三轮：去「单调梯度」表述 + 单细胞表达量口径；此前去编号"),
    (R1, "| Discussion **1,568/1,566** / Limitations **605/603** |", "| Discussion **1,589/1,587** / Limitations **605/603** |"),
    (R1, "| Conclusions + 缩略语表 | ch6 **v1.4**（删孤儿缩写 5-azaC；PPH3/IVW/dbGaP/AI 首用处补全称） |",
         "| Conclusions + 缩略语表 | ch6 **v1.6**（八层口径同步 + 删孤儿缩写 5-azaC；PPH3/IVW/dbGaP/AI 首用处补全称；第三轮：PPH0–PPH4 记法 + 表注校订） |"),
    (R1, "| Figure legends + Tables 2–7 | ch8 **v1.8**（图注去编号", "| Figure legends + Tables 2–7 | ch8 **v1.9**（图注去编号"),
    (R1, "| **9,876 含 / 9,646 不含**（= BG 762/760 + Methods 3,592/3,487 + Results 3,201/3,084 + Disc 1,568/1,566 + Lim 605/603 + Concl 148/146；2026-09-10 整改 + 补漏后统一脚本实测） |",
         "| **9,941 含 / 9,711 不含**（= BG 762/760 + Methods 3,626/3,521 + Results 3,211/3,094 + Disc 1,589/1,587 + Lim 605/603 + Concl 148/146；2026-09-10 第三轮复核后按各章增量实测） |"),
    (R1, "（本稿 **346** 达标）", "（本稿 **348** 达标）"),
    (R1, "本稿正文合计 **9,876** 含小节标题", "本稿正文合计 **9,941** 含小节标题"),
    (R1, "- v1.6（2026-09-10）：**作者四项决策全部落实**。", V23 + "- v1.6（2026-09-10）：**作者四项决策全部落实**。"),
    (R2, "*v1.7 · 2026-09-10（晚）· 阶段4.5 复现包*", R2_18 + "*v1.7 · 2026-09-10（晚）· 阶段4.5 复现包*"),
    (R2, "拼装脚本（头部装配说明 v1.7）", "拼装脚本（头部装配说明 v1.9）"),
]


def main():
    bad = 0
    for path, old, new in EDITS:
        t = open(path, encoding="utf-8").read()
        n = t.count(old)
        if n != 1:
            print(f"  FAIL hits={n}  [{os.path.basename(path)}] {old[:70]!r}")
            bad += 1
            continue
        open(path, "w", encoding="utf-8").write(t.replace(old, new, 1))
        print(f"  OK   [{os.path.basename(path)}] {old[:60]!r}")
    if bad:
        sys.exit(f"{bad} edit(s) failed")
    print("readme sync done.")


if __name__ == "__main__":
    main()
