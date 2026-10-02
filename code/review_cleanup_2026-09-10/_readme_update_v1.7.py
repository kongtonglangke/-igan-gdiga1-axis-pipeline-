#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""_readme_update_v1.7.py — 同步阶段4.5 包 README（v1.6→v1.7）与 reproducibility README（v1.1→v1.2）"""
import os
import sys

PKG = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.stdout.reconfigure(encoding="utf-8")


def patch(path, pairs):
    with open(path, encoding="utf-8", newline="") as fh:
        t = fh.read()
    for old, new, exp in pairs:
        n = t.count(old)
        if n != exp:
            raise SystemExit(f"[FAIL] {os.path.basename(path)}: {n}!= {exp}\n  >> {old[:120]}")
        t = t.replace(old, new)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write(t)
    print("OK ", os.path.basename(path), len(pairs), "pairs")


R1 = os.path.join(PKG, "README.md")
patch(R1, [
    ("# 阶段 4.5 复现包 / 投稿材料包（v1.6 · 2026-09-10）",
     "# 阶段 4.5 复现包 / 投稿材料包（v1.7 · 2026-09-10）", 1),
    ("**稿件**：Cell-type-specific genetic and epigenetic control of the galactose-deficient IgA1 axis in IgA nephropathy: a multi-omics causal framework with East Asian perspectives",
     "**稿件**：Cell-type-specific genetic and epigenetic control of the galactose-deficient IgA1 axis in IgA nephropathy: a multi-omics causal-inference framework with East Asian perspectives", 1),
    ("作者终校/建仓回填前。",
     "**本轮（v1.7）＝审稿视角全文体检后的逐条整改**（报告 v1.0 → 整改记录 v1.1）：**P0 6 项 / P1 9 项 / P2 8 项全部完成**，另修掉执行中新发现的 6 处问题（含 3 处事实性错误：Additional file 索引表 3 处图号互指错误、Abstract 整改后一度超 350 词、Fig. S1–S7 图注仍含内部编号）。三项机检通过：**正文内部编号残留 0**（仅留 Fig. 1c 图内自印的 L0–L7 行索引）、**引用首现序 [1..27] 严格递增无缺口**、**Abstract 346 ≤ 350**。作者终校/建仓回填前。", 1),
    ("    ├── code/                        M1–M6 各证据层脚本 + 主图渲染 fig_scripts + 全稿拼装脚本",
     "    ├── code/                        M1–M6 各证据层脚本 + 主图渲染 fig_scripts + 全稿拼装脚本\n    │   └── review_cleanup_2026-09-10/  审稿整改补丁脚本 5 件（可审计，非分析流程）", 1),
    ("| Title/Abstract/Keywords | ch1 v1.4（L7 入题 + 同轮瘦身） | Abstract **335**（限 350） | 全稿头部 |",
     "| Title/Abstract/Keywords | ch1 **v1.5**（审稿清洗：去编号 + 东亚句点名 rs10238682/rs7856182） | Abstract **346**（限 350 ✅） | 全稿头部 |", 1),
    ("| Background | ch2 v1.4（six→**seven** layers） | 755 含 / 754 不含 | 正文 |",
     "| Background | ch2 **v1.5**（\"heritable\" → \"influenced by common genetic variation\"） | 762 含 / 760 不含 | 正文 |", 1),
    ("| Methods（含 Table 1） | ch3 **v1.7**（L0–L7 + **East Asian regulatory cross-reference**；16 小节；S4→S4a/S4b） | **3,419 含 / 3,305 不含** | 正文 |",
     "| Methods（含 Table 1） | ch3 **v1.8**（16 小节；去编号 + 单细胞 QC 澄清 + **L4 FDR 敏感性** + 共定位先验依据 + decoupleR/IVW 口径） | **3,593 含 / 3,488 不含** | 正文 |", 1),
    ("| Results R1–R8 | ch4 **v1.6**（新增 **R8** + **R6 东亚复现句**；S4→S4a/S4b） | **3,091 含 / 2,983 不含**（8 小节） | 正文 |",
     "| Results R1–R8 | ch4 **v1.7**（去编号 + **Q1–Q4 挂标** + R1 量化 + 双分区检出率澄清；**L7 零分布改 1,000 组重跑**） | **3,203 含 / 3,086 不含**（8 小节） | 正文 |", 1),
    ("| Discussion + Limitations | ch5 **v1.5**（5 段 + **7 条**；#5 补 ImmuNexUT；长度按期刊要求**定案不瘦身**） | Discussion **1,497/1,496** / Limitations **531/530** | 正文 |",
     "| Discussion + Limitations | ch5 **v1.6**（去编号 + causal-inference 边界声明 + **三条验证清单**） | Discussion **1,568/1,566** / Limitations **605/603** | 正文 |", 1),
    ("| Conclusions + 缩略语表 | ch6 v1.3（\"Seven layers\"） | Conclusions 147 含 / 146 不含 | 正文 |",
     "| Conclusions + 缩略语表 | ch6 **v1.4**（删孤儿缩写 5-azaC；PPH3/IVW/dbGaP/AI 首用处补全称） | Conclusions **148 含 / 146 不含** | 正文 |", 1),
    ("| Title page + Declarations（8 子目） | ch7（Availability/Acknowledgements 补 L7 + **ImmuNexUT/E-GEAD-398**） | — | 投稿系统 Title page / 正文尾部 |",
     "| Title page + Declarations（8 子目） | ch7（**CRediT 重写**（并列一作自洽）+ Funding 补角色声明 + Availability 去编号） | — | 投稿系统 Title page / 正文尾部 |", 1),
    ("| Figure legends + Tables 2–7 | ch8 **v1.5**（Fig1c L7 行；**Fig2c 改成熟序图例 + panel b/a 图例诚实性修正**；Fig5 挂 S12–S13；S 子表交叉引用改写） | — | 正文尾部（图注/表） |",
     "| Figure legends + Tables 2–7 | ch8 **v1.7**（图注去编号；Fig. 1c 保留 (L0–L7) 图内索引；Fig. S6 图注同步 1,000 组零分布） | — | 正文尾部（图注/表） |", 1),
    ("| References | v1.4（**已完成一次性改号**）：**[22] ImmuNexUT / [23] CollecTRI / [24] JASPAR 2024**；原 [22]→[25] Liu、[23]→[26] Fung、[24]→[27] Kim；1–21 不动 | **27 条** | 正文尾部 |",
     "| References | **v1.5**（[27] Kim 年份 2025 → **2026**；编号体系不变） | **27 条** | 正文尾部 |", 1),
    ("| **正文合计（章节2–6，含小节标题）** | — | **9,440 含 / 9,214 不含**（= BG 755/754 + Methods 3,419/3,305 + Results 3,091/2,983 + Disc 1,497/1,496 + Lim 531/530 + Concl 147/146；原始空白分词脚本实测，取代 L7 前 6,151 的旧值） | 脚本实测定格 |",
     "| **正文合计（章节2–6，含小节标题）** | — | **9,879 含 / 9,649 不含**（= BG 762/760 + Methods 3,593/3,488 + Results 3,203/3,086 + Disc 1,568/1,566 + Lim 605/603 + Concl 148/146；2026-09-10 整改后统一脚本实测） | 脚本实测定格 |", 1),
    ("- 备注：**Discussion/Limitations 长度已定案 = 不瘦身**（2026-09-10，作者决策 ⒝「按期刊要求来」）。依据 Genome Medicine 研究论文投稿指南：唯一硬性字数上限为摘要 ≤350 词（本稿 335 达标），",
     "- 备注：**Discussion/Limitations 长度已定案 = 不瘦身**（2026-09-10，作者决策 ⒝「按期刊要求来」）。依据 Genome Medicine 研究论文投稿指南：唯一硬性字数上限为摘要 ≤350 词（本稿 **346** 达标），", 1),
    ("附：本稿正文合计 9,440 含小节标题，处于 GM 研究论文的常见区间。",
     "附：本稿正文合计 **9,879** 含小节标题；整改新增内容均为审稿人要求的透明度补强（L4 FDR 敏感性、QC 澄清、先验依据、量化句、Q1–Q4 挂标、验证清单），GM 对正文不设字数上限。", 1),
    ("v1.4（2026-09-10）：四项作者决策落地后重建——母本吸收 Fig 2c 成熟序图例、ch8 Fig 2 panel a/b 两处图例诚实性修正、S 子表交叉引用（S4a/S4b、S5a–S5d、S6a/S6b、S2a/S2b）及各章版本注；GA 段（ch1「## Graphical abstract」，不进入母本）重写为实际版面并含 L7 回路。重建后仍 **[1..27] 严格递增、无缺口、无占位符残留**。",
     "v1.4（2026-09-10）：四项作者决策落地后重建——母本吸收 Fig 2c 成熟序图例、ch8 Fig 2 panel a/b 两处图例诚实性修正、S 子表交叉引用（S4a/S4b、S5a–S5d、S6a/S6b、S2a/S2b）及各章版本注；GA 段（ch1「## Graphical abstract」，不进入母本）重写为实际版面并含 L7 回路。重建后仍 **[1..27] 严格递增、无缺口、无占位符残留**。\nv1.5（2026-09-10）：**审稿整改后重建**——全稿去内部编号、去分音符、colocalisation 英式统一、标题 causal→causal-inference；装配说明同步（去编号状态 + 补充表后缀 + 期刊字数结论）。重建后 **[1..27] 严格递增无缺口**，正文内部编号残留仅 Fig. 1c 图例的 (L0–L7)。", 1),
    ("3. **GitHub 建仓** → 上传 `reproducibility/` → Zenodo/版本 DOI → 回填 ch7 Declarations「GitHub repository DOI」与本包 README；",
     "3. **GitHub 建仓** → 上传 `reproducibility/` → Zenodo/版本 DOI → 回填 ch7 Declarations「GitHub repository DOI」与本包 README；\n\n**本轮整改新增的投稿前事项（详见整改记录第六节）**：① 署名次序/CRediT 终形态由作者确认；② 通讯作者 Yan Wang 与 [8] Wang YN 是否同一人待确认；③ GA 图内 \"upstream circuit (L7)\" 是否同去编号（改一行 + 重渲染）；④ `Additional_File_1_Figure_Legends.md` 头部 `>` 编辑版本注投稿前需删（母本的 `<!-- ASSEMBLY NOTE -->` 同理）；⑤ Funding 若投稿前立项需改资助表述。", 1),
    ("- v1.6（2026-09-10）：**作者四项决策全部落实**。",
     "- v1.7（2026-09-10）：**审稿视角全文体检后的逐条整改**。先出体检报告 `审稿视角全文检查报告_v1.0_2026-09-10.md`（只列问题、未动一字：P0 6 / P1 9 / P2 8），随后按 P0→P1→P2 全部执行，另修掉执行中新发现的 6 处问题；产出整改记录 `审稿视角全文检查报告_v1.1_整改记录_2026-09-10.md`。要点：① **全稿去内部编号**（ch3 16 个小节标题 + Results/Discussion/Declarations/图注/Additional file 共 60 余处），`(R2)` 等改图号或补充件交叉引用；② **Abstract 25→25 planned/14 computable**，并因东亚句点名变异一度超限、回修 6 词后 **346 词**；③ **单细胞 282,463 前后同值之谜**查实为 GEO 矩阵已 CellRanger 过滤、500 计数阈值未再剔除任何细胞，改为如实表述；④ **[27] Kim 年份 2025→2026**；⑤ **cover letter 六层→七层** + 补 GSE285335/ImmuNexUT；⑥ **data availability 补 4 类资源**与 ch7 对齐；⑦ naïve→naive 4 处、colocalization→colocalisation 29 处、标题 causal→causal-inference 14 处；⑧ **L4 补 BH-FDR 敏感性**（如实报出 GALNT12 罕见变异窗 q = 6.1×10⁻⁵，距全基因组显著约 3 个数量级）；⑨ **L7 经验零分布 200→1,000 重跑**（通过集合不变 36/950，emp_P 脱离分辨率下限：DNMT1 0.025 / KLF2 0.006 / SP1 0.019；motif 交集因子 TFAP2A → **MYC**，已同步 ch4 + Fig. S6）；⑩ Discussion 末段改为 **三条验证清单**；⑪ Additional file 索引表 **3 处图号互指错误**修正；⑫ 母本重建为 **v1.5**；⑬ 整改补丁脚本 5 件归档至 `reproducibility/code/review_cleanup_2026-09-10/`。\n- v1.6（2026-09-10）：**作者四项决策全部落实**。", 1),
])

R2 = os.path.join(PKG, "reproducibility", "README.md")
patch(R2, [
    ("| `code/_build_submission_fulltext.py` | 章节定稿 → 单文件全稿（Manuscript）拼装脚本 | 全稿组装 |",
     "| `code/_build_submission_fulltext.py` | 章节定稿 → 单文件全稿（Manuscript）拼装脚本（v1.5 头部同步去编号状态） | 全稿组装 |", 1),
    ("| `code/_tmp_geo_check.py` | GEO 系列元数据/文献核对 | Methods 数据决议 |",
     "| `code/_tmp_geo_check.py` | GEO 系列元数据/文献核对 | Methods 数据决议 |\n| `code/sc_regulatory/` | 单细胞 atlas（step1）、**上游调控子推断（step2；`--k-null` 默认 1,000）**、甲基化敏感 motif（step3） | Results R8；Fig. S5–S6；Table S9/S10 |\n| `code/directional_validation/` | GSE285335 载入与 B 细胞亚群表达/轨迹检验 | Results R2 单细胞交叉核对；Fig. S4；Table S7/S8 |\n| `code/C1_eas_eqtl/` | ImmuNexUT E-GEAD-398 东亚调控交叉参照（step1–3） | Results R6；Table S12/S13 |\n| `code/review_cleanup_2026-09-10/` | **审稿整改补丁脚本 5 件**（体例统一、去内部编号、定点改写、L7 数字同步、Additional file 图注修正）；仅供审计，不参与分析流程 | — |", 1),
    ("- `00_rawdata` 原样下载缓存**不随本包分发**（体积/许可），由下载脚本与 accession 重建。",
     "- `00_rawdata` 原样下载缓存**不随本包分发**（体积/许可），由下载脚本与 accession 重建。\n- 单细胞调控子推断的环境：`C:/Users/user/.workbuddy/binaries/python/envs/scRNA`（numpy 2.5.3 / scipy / pandas / statsmodels 0.15 / scikit-learn）。\n- **2026-09-10 重跑记录**：`step2_regulon.py --k-null 1000`（原 200）——经验零分布由 200 组目标数匹配随机基因集提升到 **1,000 组**，使经验 P 脱离 1/K 分辨率下限；**通过集合不变（36 / 950，3.8%）**，DNMT1 0.010→0.025、KLF2 0.005→0.006、SP1 0.015→0.019；下游 step3 交集因子随之由 TFAP2A 变为 **MYC**（ch4 R8 与 Fig. S6 已同步）。脚本新增 `--k-null` 参数（默认 1,000），`step2_summary.json` 记录 `k_null`。", 1),
    ("全表索引（S1–S13）见 `../submission/additional_file_1/Additional_File_1_Figure_Legends.md`。本 `data/` 目录保留**分析内部命名**（不含 Sxx 后缀），两者为\"分析源 → 投稿件\"映射关系，不做重命名以免断脚本。",
     "全表索引（S1–S13）见 `../submission/additional_file_1/Additional_File_1_Figure_Legends.md`。本 `data/` 目录保留**分析内部命名**（不含 Sxx 后缀），两者为\"分析源 → 投稿件\"映射关系，不做重命名以免断脚本。\n\n> **关于本文件与正文的编号口径差异**：上表\"对应稿件\"列沿用分析期的内部证据层编号（L0–L7 / R1–R8）作为**代码与结果的定位索引**；2026-09-10 审稿整改后，**投稿正文已不再使用这些编号**（改为图号/补充件交叉引用），本文件是这些编号唯一保留的地方。", 1),
    ("*v1.1 · 2026-09-10 · 阶段4.5 复现包*",
     "*v1.2 · 2026-09-10 · 阶段4.5 复现包* — v1.1→v1.2：补列 `sc_regulatory/`、`directional_validation/`、`C1_eas_eqtl/`、`review_cleanup_2026-09-10/` 四个子目录；登记 step2 `--k-null 1000` 重跑与 motif 交集因子变化；说明内部编号仅存于本文件的定位口径。", 1),
])
print("done")
