# -*- coding: utf-8 -*-
"""
M6-1 干预映射证据输入表构建
============================
把阶段2评估的 5 个 O-糖基化轴基因（C1GALT1 / C1GALT1C1 / GALNT2 / GALNT12 / ST6GALNAC2）
映射到表观遗传药物 / 天然产物（DNMTi / HDACi / EGCG），产出：
  1) M6_药物目录.tsv            —— 候选干预药物/分子目录（含实际表观靶酶、批准状态、B/浆细胞证据）
  2) M6_靶点药物证据表.tsv      —— 5 轴基因 × {功能/表达状态/表观窗口/我方遗传证据/干预方向/候选药/证据等级/直接可成药性/评分}
证据来源口径（详见 M6_方法学说明段落_20260903.md）：
  - EpiFactors：表观酶角色分类（DNMT1=writer; HDAC1/2/3/6=eraser）——按数据库公开分类标注
  - DrugBank：批准状态按公开记录（FDA/NMPA 批准年份与适应症）
  - Pharos (IDG)：直接可成药性等级 Tclin/Tchem/Tbio/Tdark（2026-09-03 实测）
  - CMap：连接性查询未执行（许可限制），表中药物-基因关系基于文献 + 公开注释，已在说明中声明
本脚本为「纯编制 + 可复现」结构：所有内容以 Python 字面量沉淀，改动即重跑即可。
"""
import os
import pandas as pd

OUT = os.path.dirname(os.path.abspath(__file__))

# =====================================================================
# 1. 药物目录（10 条：天然产物 + DNMTi + HDACi；另附范围外临床对照注释行）
# =====================================================================
DRUGS = [
    dict(
        agent="EGCG（表没食子儿茶素没食子酸酯）",
        class_="天然多酚（绿茶）",
        target_enzyme="DNMT1/DNMT3B（直接结合抑制，Ki≈6.9 µM，非核苷类）；HDAC1/2（可逆抑制）；对 HAT 亦有报道",
        epifactors_role="上游抑制 writer(DNMT)/eraser(HDAC)，属天然表观调节剂",
        modality="非核苷类；非细胞毒机制",
        approved="未作为处方药上市（外用制剂 Veregen/Polyphenon E 用于疣；膳食/研究级）",
        beell_plasma_evidence="对 B 淋巴细胞具免疫调节作用（综述级）；体外可再激活甲基化沉默基因；论文1 EGCG 细胞/动物功能结果（引用待预印本链接回填）",
        usage_hypothesis="候选『转录恢复』天然分子：假设上调 B/浆细胞 C1GALT1/C1GALT1C1 转录→降低 Gd-IgA1；安全性最友好的临床前手柄",
        source_note="Fang 2003（Ki 6.89 µM）；Negri 2018 Nutrients 综述；HeLa 直接结合 DNMT1/DNMT3B/HDAC1；现代免疫学 2018 综述",
        risk_note="浓度依赖：高浓度促氧化/细胞毒；口服生物利用度低；多靶点混杂",
    ),
    dict(
        agent="Azacitidine（5-氮杂胞苷，5-aza-CR；Vidaza）",
        class_="核苷类 DNMT 抑制剂",
        target_enzyme="DNMT1（掺入 RNA/DNA 后被 DNMT 捕获→降解→被动去甲基化）",
        epifactors_role="writer 抑制剂",
        modality="核苷类（细胞毒背景）",
        approved="FDA 2004 批准（MDS 等髓系肿瘤）",
        beell_plasma_evidence="体外可逆转 IL-4/IL-17 诱导的 IgA1 低半乳糖基化（IgAN B 细胞/细胞系），恢复 C1GALT1C1 表达——IgAN 背景直接功能证据（引用见 T 细胞综述）",
        usage_hypothesis="机制探针/概念验证药：验证『去甲基化恢复 C1GALT1C1/C1GALT1→降低 Gd-IgA1』因果链",
        source_note="Clin Exp Nephrol 2019（T cells in IgAN 综述）及其引用的 IL-4/IL-17 原始研究",
        risk_note="全身骨髓抑制/肝毒性/致畸，不适直接入药；仅作窗口机制验证",
    ),
    dict(
        agent="Decitabine（5-氮杂-2'-脱氧胞苷，5-aza-CdR；Dacogen）",
        class_="核苷类 DNMT 抑制剂",
        target_enzyme="DNMT1（DNA 掺入，捕获并耗竭 DNMT）",
        epifactors_role="writer 抑制剂",
        modality="核苷类（细胞毒背景）",
        approved="FDA 2006 批准（MDS/AML）",
        beell_plasma_evidence="同 5-aza-CR（同一机制家族）；IgAN B 细胞直接证据弱于 5-aza-CR 文献",
        usage_hypothesis="同 5-aza-CR，作为去甲基化窗口的替代工具药",
        source_note="公开药理学记录；IgAN 应用为外推假设",
        risk_note="同核苷类毒性；长期去甲基化存在二次肿瘤风险",
    ),
    dict(
        agent="RG108",
        class_="非核苷类 DNMT1 抑制剂（研究用）",
        target_enzyme="DNMT1（活性位点直接结合，非掺入型）",
        epifactors_role="writer 抑制剂",
        modality="非核苷类",
        approved="未批准（实验室工具药）",
        beell_plasma_evidence="无 IgAN/B 细胞直接证据",
        usage_hypothesis="低毒去甲基化探针备选（若需规避核苷类毒性）",
        source_note="公开药理学记录",
        risk_note="体外活性弱、细胞摄取有限",
    ),
    dict(
        agent="Vorinostat（伏立诺他，SAHA；Zolinza）",
        class_="羟肟酸类 HDAC 抑制剂（class I/II 广谱）",
        target_enzyme="HDAC1/2/3/6 等（Zn2+ 依赖）",
        epifactors_role="eraser 抑制剂",
        modality="小分子（广谱表观）",
        approved="FDA 2006 批准（皮肤 T 细胞淋巴瘤 CTCL）",
        beell_plasma_evidence="小鼠：损害初次抗体应答但保留记忆 B 细胞（Waibel 2015 Nat Commun）——支持 B 细胞应答窗口可调",
        usage_hypothesis="若 C1GALT1C1/C1GALT1 沉默含去乙酰化成分，class I HDACi 或可再激活（假设，需功能检验）；亦代表『浆细胞/抗体输出调节』方向",
        source_note="Waibel et al. Nat Commun 6:7838 (2015)",
        risk_note="骨髓抑制/腹泻/心脏 QT 风险；广谱表观脱靶",
    ),
    dict(
        agent="Romidepsin（罗米地辛；Istodax）",
        class_="环肽类 HDAC 抑制剂（class I 选择性：HDAC1/2）",
        target_enzyme="HDAC1/HDAC2（class I）",
        epifactors_role="eraser 抑制剂",
        modality="小分子（class I 选择性）",
        approved="FDA 2009 CTCL / 2011 PTCL",
        beell_plasma_evidence="B/浆细胞肿瘤敏感；IgAN 无直接证据",
        usage_hypothesis="class I（HDAC1/2）选择性——若轴基因表观窗口由 HDAC1/2 介导则选择性更优（假设）",
        source_note="公开药理学记录；MM 综述（HDACi 在浆细胞肿瘤）",
        risk_note="QT/骨髓毒性；静脉给药",
    ),
    dict(
        agent="Panobinostat（帕比司他，LBH589；Farydak）",
        class_="羟肟酸类泛-HDAC 抑制剂（class I/II/IV，class I 效力更优）",
        target_enzyme="HDAC1/2/3/6 等",
        epifactors_role="eraser 抑制剂",
        modality="小分子（广谱）",
        approved="FDA 2015 批准（多发性骨髓瘤，与硼替佐米+地塞米松联用）",
        beell_plasma_evidence="显著削减自身反应性浆细胞、自身抗体与肾炎（MRL/lpr 狼疮小鼠，Waibel 2015）；多发性骨髓瘤（浆细胞肿瘤）临床获批背景",
        usage_hypothesis="『浆细胞输出/自身抗体削减』方向的动物概念证据最强者；理论上可压低 Gd-IgA1 产生源（非转录恢复方向）",
        source_note="Waibel et al. Nat Commun 6:7838 (2015)；FDA 2015",
        risk_note="显著骨髓抑制/消化道毒性；广谱脱靶；不提示直接用于 IgAN",
    ),
    dict(
        agent="Entinostat（恩替诺特，MS-275）",
        class_="苯酰胺类 HDAC 抑制剂（class I 选择性：HDAC1/2/3）",
        target_enzyme="HDAC1/2/3（class I）",
        epifactors_role="eraser 抑制剂",
        modality="小分子（class I 选择性，口服）",
        approved="未批准（曾进实体瘤/血液肿瘤多期临床）",
        beell_plasma_evidence="无 IgAN 直接证据；class I 谱与 C1GALT1C1 表观窗口假设较契合",
        usage_hypothesis="class I 选择性口服探针备选（若获批前景存在）",
        source_note="公开药理学记录",
        risk_note="临床开发中安全性数据不完整",
    ),
    dict(
        agent="Chidamide（西达本胺；Epidaza）",
        class_="苯酰胺类 HDAC 抑制剂（class I 选择性 HDAC1/2/3/10）",
        target_enzyme="HDAC1/2/3/10（class I）",
        epifactors_role="eraser 抑制剂",
        modality="小分子（class I 选择性，口服）",
        approved="中国 NMPA 2015 批准（外周 T 细胞淋巴瘤 PTCL）；表观药物中国本土可及代表",
        beell_plasma_evidence="无 IgAN 直接证据",
        usage_hypothesis="中国语境下『class I HDACi 可及性』的代表；假设同上",
        source_note="NMPA 公开记录",
        risk_note="血液学毒性；超适应症使用需伦理/法规评估",
    ),
    dict(
        agent="Belinostat（贝利司他；Beleodaq）",
        class_="羟肟酸类 HDAC 抑制剂（class I/II）",
        target_enzyme="HDAC1/2/3/6 等",
        epifactors_role="eraser 抑制剂",
        modality="小分子（广谱）",
        approved="FDA 2014 批准（PTCL）",
        beell_plasma_evidence="无 IgAN 直接证据",
        usage_hypothesis="同类补充（静脉），逻辑同 vorinostat",
        source_note="公开药理学记录",
        risk_note="同羟肟酸类毒性谱",
    ),
    dict(
        agent="Trichostatin A（曲古抑菌素 A，TSA）",
        class_="羟肟酸类泛-HDAC 抑制剂（工具药）",
        target_enzyme="HDAC1/2/3/6 等（class I/II）",
        epifactors_role="eraser 抑制剂",
        modality="小分子（研究工具）",
        approved="未批准（经典体外工具药）",
        beell_plasma_evidence="HDACi 领域机制研究参照（体外）",
        usage_hypothesis="体外『HDAC 窗口』探针",
        source_note="公开记录",
        risk_note="体内不稳定；仅体外",
    ),
]

DRUG_COLS = ["agent","class_","target_enzyme","epifactors_role","modality","approved",
             "beell_plasma_evidence","usage_hypothesis","source_note","risk_note"]
DRUG_RENAME = {"agent":"药物/分子","class_":"类别","target_enzyme":"实际靶酶(表观)",
               "epifactors_role":"EpiFactors 角色","modality":"化学型/毒性特征",
               "approved":"批准状态(公开记录)","beell_plasma_evidence":"B/浆细胞相关证据",
               "usage_hypothesis":"在 IgAN 轴上的潜在用途(假设)","source_note":"证据出处",
               "risk_note":"风险/备注"}

# =====================================================================
# 2. 5 轴基因靶点-药物证据行
# =====================================================================
GENES = [
    dict(
        gene="C1GALT1", chr="7",
        node="IgA1 铰链区 O-半乳糖基化核心酶（Core1 β1,3-半乳糖基转移酶，T-synthase）；生成 Galβ1-3GalNAc（T 抗原）",
        igan_expr_state="IgAN PBMC 显著下调（GSE73953 logFC=-1.73, P=6.3e-11，bulk 混合群、n_ctrl=2 谨慎）；肾小球亦下调（GSE93798 P=3.4e-5）；文献：IgAN 的 Gd-IgA1 生成细胞 C1GALT1 表达/活性降低",
        window_acquired="获得性层：IL-4/IL-17/TGF-β 下调 C1GALT1 mRNA（机制多认为经转录/表观，但直接启动子甲基化因果证据主要落在伴侣基因 C1GALT1C1）；跨癌种 TCGA 甲基化-表达负相关见于肾/乳腺/前列腺等癌（相关性而非因果，且癌型依赖——GI 癌反而低甲基化高表达）",
        window_germline="种系层（我方 M4+阶段2）：lead rs13226913/rs10238682 是 C1GALT1 基因体内/3' 侧 5 个 CpG 的强 cis-mQTL（P≤1e-321，基因体甲基化与转录同向耦合，非启动子窗）；已发表 QT-GWAS 提示核心启动子 SNP rs758263 与 miRNA(rs1047763/miR-148b) 调节",
        anchor_ours="OneK1K B_naive eQTL P=9.4e-6（复现）；GTEx blood 5.6e-12 / kidney 1.9e-5；coloc×IgAN 全阴性（max PPH4=0.011）；M4 mQTL×IgAN coloc PPH4≤0.014（甲基化侧亦不共享）",
        anchor_disease="IgAN(GWAS)：rs13226913 P=0.284、rs10238682 P=0.222（阴性）；Gd-IgA1(Wang2021 表型)：rs10238682 P=1.2e-9、rs13226913 P=0.011（机制表型强锚）",
        goal="恢复/上调 B/浆细胞 C1GALT1 转录→提高 IgA1 铰链区半乳糖基化→降低 Gd-IgA1",
        candidates="DNMTi(5-aza-CR/地西他滨/RG108) 与 EGCG（非核苷低毒）：假设靶向『转录沉默窗』（跨疾病甲基化-表达负相关为相关性证据，无 C1GALT1 特异性去甲基化功能因果）；class I HDACi(vorinostat/entinostat/西达本胺)：假设针对去乙酰化沉默窗；论文1 EGCG 功能衔接（引用待回填）——注意多数肿瘤中 C1GALT1 反而上调（促癌），上调干预须 B 细胞/组织特异",
        ev_tier="C 级（假设为主）：IgAN B 细胞『药物恢复 C1GALT1』直接功能证据缺；跨癌种甲基化-表达数据为相关性（非因果）",
        pharos="Tbio（0 配体/0 批准药）→ 直接酶位点不可成药；唯一现成手柄=转录/表观药物",
        risk="全身性 DNMTi/HDACi 毒性；表观广谱脱靶；多数实体瘤中 C1GALT1 过表达为促癌/不良预后面（上调干预需 B 细胞向精准，避免系统性促癌风险）",
        G=3.0, D=1.5, F=1.0, P=1.0,
        g_just="B 系(naive)eQTL 复现(9.4e-6) + Gd-IgA1 表型强锚(1.2e-9) + 表达方向与机制一致（基因体内 cis-mQTL 亦强阳性）",
        d_just="转录窗口明确（M4 基因体 mQTL×eQTL 同向耦合 + 跨癌种甲基化-表达负相关为相关性），DNMTi/EGCG 手柄理论可及；但 IgAN B 细胞直接恢复证据缺、C1GALT1 特异性启动子甲基化因果未确立 → 低于 C1GALT1C1 的确定性",
        f_just="跨疾病(肿瘤)甲基化-表达相关性（非功能因果）；无 IgAN 功能数据（论文1 EGCG 结果未回填暂不计分）",
        p_just="核苷类 DNMTi 全身毒性 + 广谱表观脱靶 + 肿瘤促癌矛盾面",
    ),
    dict(
        gene="C1GALT1C1", chr="X",
        node="C1GALT1 特异性分子伴侣（COSMC，T-synthase 伴侣）；保证 C1GALT1 折叠/稳定/活性",
        igan_expr_state="IgAN PBMC 极强下调（GSE73953 logFC=-2.46, P=8.5e-10）；文献：Gd-IgA1 生成细胞 C1GALT1C1 表达降低；IL-4 对 IgAN B 细胞的 C1GALT1C1 下调更强（IgAN B 细胞对 IL-4 更敏感）",
        window_acquired="获得性层（最强证据）：IL-4/IL-17 诱导 C1GALT1C1 启动子 CpG 岛超甲基化→mRNA 下调→分泌异常糖基化 IgA1；5-氮杂胞苷可逆转（IgAN B 细胞/细胞系体外功能证据）——『可逆的获得性表观沉默窗』成立",
        window_germline="种系层：位于 X 染色体，cis-eQTL 物理不可测（阶段2 诚实缺口）；Wang2021 Gd-IgA1 关联 lead rs5910940（X:120.6Mb, P=0.004）",
        anchor_ours="阶段2：X 染色体 cis-eQTL 不可测/未命中；coloc 不可算；M5：rs5910940 在东亚 Gd-IgA1 P=0.004（无 IgAN 对照数据）——遗传锚点弱于 C1GALT1，生物学锚点强",
        anchor_disease="IgAN：不可测（X）；Gd-IgA1：rs5910940 P=0.004（名义）",
        goal="去甲基化恢复 C1GALT1C1 表达→恢复 C1GALT1 活性→降低 Gd-IgA1",
        candidates="DNMTi 直接命中：5-氮杂胞苷（IgAN B 细胞功能证据）、地西他滨、非核苷 EGCG/RG108（低毒备选）；HDACi 证据弱",
        ev_tier="A 级（5-aza-CR 在 IgAN B 细胞/细胞系逆转 IL-4/IL-17 低糖基化，体外功能直接证据，经 T 细胞综述转引原始文献）；体内证据缺",
        pharos="Tbio（0 配体/0 批准药）→ 直接不可成药；药物作用对象为上游 DNMT1（writer）而非伴侣本身",
        risk="核苷类 DNMTi 全身毒性（骨髓/肝/致畸）；启动子去甲基化广谱性；X 染色体剂量/失活混杂",
        G=2.0, D=3.0, F=2.0, P=1.0,
        g_just="X 染色体 cis-eQTL 不可测，遗传锚点弱；Gd-IgA1 名义关联(0.004) + IgAN B 细胞表达/IL-4 敏感性文献构成机制锚（遗传分项下调）",
        d_just="5-aza-CR 在 IgAN B 细胞直接逆转 C1GALT1C1 启动子超甲基化沉默——全部 5 基因中唯一『IgAN 背景直接功能证据』的药物-窗口匹配",
        f_just="IgAN 来源 B 细胞/细胞系体外功能证据（5-aza-CR 逆转）",
        p_just="核苷类 DNMTi 全身毒性 + 去甲基化广谱 + 无体内验证",
    ),
    dict(
        gene="GALNT12", chr="9",
        node="多肽 GalNAc 转移酶（GalNAc-T12）：铰链区 O-GalNAc 起始酶（初始糖基化步骤成员）",
        igan_expr_state="IgAN PBMC 显著下调（GSE73953 logFC=-2.12, P=5.6e-4，bulk 谨慎）；肾组织无显著差异",
        window_acquired="无 IgAN 特异表观窗口证据",
        window_germline="种系层：Wang2021 Gd-IgA1 lead rs7856182（chr9:98.87Mb，距基因约 19kb，非基因体内）关联 P=2.4e-9（机制表型强锚），但 B 系/全血 cis-eQTL 均不显著（阶段2 诚实缺口，转 V1 验证关注项）",
        anchor_ours="OneK1K B 系/GTEx blood/kidney eQTL 均 P>0.1（B_naive 0.246/记忆 0.36 等不显著）；coloc×IgAN 阴性（PPH4 max 0.0051）；cis-eQTL 区域顶变异反查 IgAN 未达基因组显著",
        anchor_disease="IgAN：rs7856182 P=0.100（阴性）；Gd-IgA1：rs7856182 P=2.4e-9（机制表型强锚）",
        goal="（若 V1 验证支持）恢复 B 细胞 GALNT12 表达以保初始 GalNAc 添加；目前证据不足以支撑干预主张",
        candidates="无直接药物（Pharos Tdark，0 化合物）；表观药物无窗口依据——本基因当前定位为『验证包关注项』而非干预锚点",
        ev_tier="C 级/证据不足：机制表型(Gd-IgA1)遗传强锚，但 B 系 eQTL 阴性 → 细胞级转录窗口未确立，干预主张不成立",
        pharos="Tdark（0 配体/0 批准药）→ 不可成药",
        risk="即使未来验证阳性，GalNAc-T 家族(20 同工酶)选择性抑制/激活均困难（同源结构相似）",
        G=2.0, D=0.5, F=0.0, P=0.0,
        g_just="Gd-IgA1 机制表型强锚(2.4e-9)但 B 系/全血 eQTL 阴性——遗传锚为『表型→基因』外推，细胞级证据缺，故给 2.0 上限内折中",
        d_just="Pharos Tdark 无化合物；无表观窗口证据（仅象征性 0.5 表示理论可及转录水平干预）",
        f_just="无任何功能药理学证据",
        p_just="暂无具体干预管线，无对应风险项",
    ),
    dict(
        gene="GALNT2", chr="1",
        node="多肽 GalNAc 转移酶（GalNAc-T2）：铰链区 O-GalNAc 起始酶；UniProt 标注『可能参与 IgA1 铰链区 O-糖基化』；亦调控 ANGPTL3/ApoC-III(脂代谢)与胰岛素受体",
        igan_expr_state="IgAN PBMC 显著下调（GSE73953 logFC=-4.43, P=4.6e-6，最强下调，bulk 谨慎）；肾组织反而上调（GSE115857/93798）——组织方向分离",
        window_acquired="无 IgAN 特异表观窗口证据",
        window_germline="种系层：无 Gd-IgA1/IgAN 直接 lead；区域 eQTL 顶变异反查 IgAN 未显著；M2b SuSiE 复核 GALNT2×B_naive 有官方 credible set 但逐信号 PPH4≤0.031（阴性）",
        anchor_ours="coloc×IgAN 5 语境 max PPH4=0.0090（不支持共享）；表达锚点主要来自 M1 GEO（bulk）",
        anchor_disease="IgAN/Gd-IgA1：无直接 lead 关联记录（较 C1GALT1/GALNT12 弱）",
        goal="证据不足以支撑 IgAN 干预主张；如未来证明其参与 IgA1 铰链糖基化且方向明确，可考虑 GalNAc-T2 选择性调节（研究级配体）",
        candidates="GalNAc-T2 研究级配体（Pharos Tchem，3 活性配体；最新双域抑制剂 IC50≈21 µM，选择性优于 T1/T3——JACS Au 2024）；表观药物无窗口依据；当前定位『表达锚点基因』非干预锚点",
        ev_tier="C 级/证据不足：脂代谢与 CDG 病证据为主，IgAN 方向无遗传 lead、无功能药理证据",
        pharos="Tchem（3 研究配体，0 批准药）→ 无临床直接药",
        risk="GalNAc-T2 缺失=GALNT2-CDG（神经发育病）；全身抑制有害；脂代谢旁路效应（ANGPTL3）",
        G=1.5, D=1.0, F=0.0, P=0.0,
        g_just="M1 GEO 多数据集表达方向(强下调)+UniProt 铰链糖基化注释构成机制线索，但无 Gd-IgA1/IgAN 遗传 lead、B 系 eQTL 证据弱",
        d_just="仅研究级 GalNAc-T2 配体（非表观药物窗口）；表观手柄无依据",
        f_just="无功能药理学证据",
        p_just="暂无具体干预管线；若推进需防 GALNT2-CDG 表型/脂代谢脱靶",
    ),
    dict(
        gene="ST6GALNAC2", chr="17",
        node="α-2,6-唾液酸转移酶（ST6GalNAc-II）：对铰链区 GalNAc 过早唾液酸化→阻止半乳糖添加→促进 Gd-IgA1 生成（『过早唾液酸化』窗口）",
        igan_expr_state="IgAN PBMC 显著下调（GSE73953 logFC=-2.16, P=5.1e-10）；文献：Gd-IgA1 生成细胞 ST6GALNAC2 表达升高（与 C1GALT1 降低并存）；肾活检上调（GSE115857）——需注意 bulk 组织方向分离与文献『B 细胞中升高』不完全一致（亚群特异）",
        window_acquired="获得性层：IL-6 上调 ST6GALNAC2 活性/表达（IgA1 分泌细胞）；无甲基化/乙酰化特异证据",
        window_germline="种系层：无 Gd-IgA1/IgAN 直接 lead；M2b SuSiE：ST6GALNAC2×blood 官方 credible set L1 逐信号 PPH4≈0.031（阴性）；区域 eQTL 顶变异反查 IgAN 未显著",
        anchor_ours="coloc×IgAN max PPH4=0.0308（不支持共享）；表达/活性方向由文献（IL-6→上调）与 M1 GEO（PBMC 下调，方向矛盾提示亚群特异）支持",
        anchor_disease="IgAN/Gd-IgA1：无直接 lead 关联记录",
        goal="（假设）抑制 ST6GALNAC2 过早唾液酸化→保留 GalNAc 供半乳糖基化→降低 Gd-IgA1；注意需在 B 细胞亚群内先验证『升高』方向（与 PBMC bulk 下调矛盾）",
        candidates="唾液酸转移酶直接抑制剂仅研究级（CMP-Neu5Ac 过渡态类似物：DANA 类，细胞渗透/磷酸酶稳定性差，无临床）；无批准药（Pharos Tbio）；表观药物（DNMTi/HDACi）调控无依据",
        ev_tier="C 级/假设：唾液酸转移酶抑制概念来自肿瘤糖生物学，IgAN 无功能药理证据；方向待单细胞亚群确认",
        pharos="Tbio（0 配体/0 批准药）→ 不可成药",
        risk="唾液酸化广谱抑制影响免疫细胞迁移/黏附；方向(上调 or 下调)存在 bulk vs 文献矛盾，需先解",
        G=1.5, D=1.0, F=0.0, P=0.0,
        g_just="文献（IL-6→上调，促 Gd-IgA1）+M1 PBMC 显著下调（方向矛盾提示亚群特异）构成机制线索；无遗传 lead",
        d_just="直接抑制概念存在（研究级唾液酸转移酶抑制剂）但不可成药性高；表观手柄无依据",
        f_just="无功能药理学证据",
        p_just="唾液酸化广谱脱靶风险（免疫细胞功能）；方向未定",
    ),
]

GENE_COLS = ["gene","chr","node","igan_expr_state","window_acquired","window_germline",
             "anchor_ours","anchor_disease","goal","candidates","ev_tier","pharos","risk",
             "G","D","F","P","interv","g_just","d_just","f_just","p_just"]
GENE_RENAME = {"gene":"基因","chr":"染色体","node":"蛋白功能/通路节点",
               "igan_expr_state":"IgAN 中表达状态(文献+M1 GEO)",
               "window_acquired":"获得性表观窗口(可干预锚点)",
               "window_germline":"种系调控窗口(我方遗传证据)",
               "anchor_ours":"我方遗传/表观证据(阶段1-3)",
               "anchor_disease":"疾病/表型锚(Gd-IgA1 vs IgAN)",
               "goal":"期望干预方向",
               "candidates":"候选干预映射",
               "ev_tier":"功能证据等级(A/B/C)",
               "pharos":"直接可成药性(Pharos IDG, 2026-09-03)",
               "risk":"特异/安全风险",
               "G":"评分_遗传锚点G(0-3)","D":"评分_药物手柄D(0-3)",
               "F":"评分_功能证据F(0-2)","P":"评分_风险扣分P(0-2)","interv":"可干预性总分(0-8)",
               "g_just":"G 评分依据","d_just":"D 评分依据","f_just":"F 评分依据","p_just":"P 评分依据"}

def score_row(r):
    r = dict(r)
    r["interv"] = round(r["G"] + r["D"] + r["F"] - r["P"], 1)
    return r

def main():
    df_drug = pd.DataFrame([{k: d[k] for k in DRUG_COLS} for d in DRUGS])
    df_drug = df_drug.rename(columns=DRUG_RENAME)
    df_drug.to_csv(os.path.join(OUT, "M6_药物目录.tsv"), sep="\t", index=False)

    rows = [score_row(g) for g in GENES]
    df_g = pd.DataFrame([{k: r[k] for k in GENE_COLS} for r in rows])
    df_g = df_g.rename(columns=GENE_RENAME)
    df_g.to_csv(os.path.join(OUT, "M6_靶点药物证据表.tsv"), sep="\t", index=False)

    # 评分简表（供出图与快速核对）
    df_score = df_g[[ "基因", "评分_遗传锚点G(0-3)", "评分_药物手柄D(0-3)",
                      "评分_功能证据F(0-2)", "评分_风险扣分P(0-2)", "可干预性总分(0-8)"]].copy()
    df_score = df_score.sort_values("可干预性总分(0-8)", ascending=False).reset_index(drop=True)
    df_score.to_csv(os.path.join(OUT, "M6_可干预性评分.tsv"), sep="\t", index=False)

    print("药物目录:", df_drug.shape)
    print("靶点-药物证据表:", df_g.shape)
    print(df_score.to_string(index=False))
    print("OK ->", OUT)

if __name__ == "__main__":
    main()
