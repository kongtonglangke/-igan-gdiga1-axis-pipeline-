# -*- coding: utf-8 -*-
"""阶段2 主证据1 —— 步骤B v2：lead Gd-IgA1 位点 × 细胞型 eQTL 效应矩阵构建
输入：步骤A 落盘的 region_cache/{qtd}.{chrom}.tsv（all.tsv.gz 区域原始行，19 列）
输出：
  1) lead_hits_all.tsv    —— 4 lead SNP 本尊在所有语境中的精确位置命中
  2) lead_eQTL_matrix.tsv —— 汇总矩阵（修复：float 比较取 min-p 行，附 beta/se）
  3) lead_region_top.tsv  —— 各 lead 区域最强 cis-eQTL（上下文）
列结构（已实测，19 列）：
  0 molecular_trait_id 1 chromosome 2 position 3 ref 4 alt 5 variant 6 ma_samples
  7 maf 8 pvalue 9 beta 10 se 11 type 12 ac 13 an 14 r2 15 molecular_trait_object_id
  16 gene_id 17 median_tpm 18 rsid
"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import os, glob

CACHE = _paths.at("00_rawdata/eQTL_Catalogue/region_cache")
OUT = _paths.at("阶段2/主证据1_lead_eQTL")
os.makedirs(OUT, exist_ok=True)

QTD2LABEL = {
    "QTD000608": "OneK1K_B_naive",
    "QTD000607": "OneK1K_B_memory",
    "QTD000606": "OneK1K_B_intermediate",
    "QTD000356": "GTEx_v10_blood",
    "QTD000261": "GTEx_v10_kidney_cortex",
}

LEADS = [
    ("rs13226913", "7", 7207215,   "C1GALT1",  "ENSG00000106392"),
    ("rs10238682", "7", 7215386,   "C1GALT1",  "ENSG00000106392"),
    ("rs7856182",  "9", 98871030,  "GALNT12",  "ENSG00000119514"),
    ("rs5910940",  "X", 120624955, "C1GALT1C1","ENSG00000171155"),
]
# mygene 实测 GRCh38 基因体坐标（2026-09-02）
GENE_BODY = {
    "C1GALT1":   ("7", 7156934, 7259333),
    "C1GALT1C1": ("X", 120614846, 120630115),
    "GALNT12":   ("9", 98807613, 98851267),
    "GALNT2":    ("1", 230057990, 230282128),
    "ST6GALNAC2":("17", 76564424, 76586956),
}
AXIS = {v: k for k, v in {
    "C1GALT1": "ENSG00000106392", "C1GALT1C1": "ENSG00000171155",
    "GALNT12": "ENSG00000119514", "GALNT2": "ENSG00000143641",
    "ST6GALNAC2": "ENSG00000070731"}.items()}

import json
_SYM = {}
_sym_path = _paths.legacy(".workbuddy/tmp/tabix/ensg2sym.json")
if os.path.exists(_sym_path):
    with open(_sym_path, encoding="utf-8") as f:
        _SYM = json.load(f)
def sym(gid):
    return _SYM.get(gid, AXIS.get(gid, ""))

def parse_rows(path):
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for ln in f:
            ln = ln.rstrip("\r\n")
            if not ln.strip() or ln.startswith("molecular_trait_id"):
                continue
            yield ln.split("\t")

def main():
    all_hits, region_top = [], []
    for path in sorted(glob.glob(CACHE("*.tsv"))):
        qtd = os.path.basename(path).split(".")[0]
        chrom = os.path.basename(path).split(".")[1]
        label = QTD2LABEL.get(qtd, qtd)
        leads_here = [L for L in LEADS if L[1] == chrom]
        rows = [f for f in parse_rows(path) if len(f) >= 19]
        # lead 本尊精确位置命中
        for rsid, lchrom, lpos, gsymbol, gensg in leads_here:
            hit = [f for f in rows if f[2] == str(lpos)]
            if not hit:
                all_hits.append([rsid, lchrom, str(lpos), gsymbol, label, qtd,
                                 "NO_LEAD_SNP_IN_REGION", "", "", "", "", "", "", "", ""])
                continue
            for f in hit:
                all_hits.append([rsid, lchrom, str(lpos), gsymbol, label, qtd,
                                 f[0], AXIS.get(f[0], ""), f[8], f[9], f[10], f[7], f[18], f[5], f[14]])
        # 区域 top（每区域前 25 最显著）
        tmp = [(float(f[8]), f) for f in rows]
        tmp.sort(key=lambda x: x[0])
        for pv, f in tmp[:25]:
            region_top.append([chrom, label, qtd, f[0], sym(f[0]),
                               f[2], f[8], f[9], f[10], f[18], f[5], f[14]])

    with open(OUT("lead_hits_all.tsv"), "w", encoding="utf-8") as o:
        o.write("lead_rsid\tlead_chr\tlead_pos\tlead_gene_symbol\tdataset\tdataset_qtd\t"
                "eQTL_gene_id\teQTL_gene_symbol\tpvalue\tbeta\tse\tmaf\trsid_col\tvariant\tLD_r2\n")
        for r in all_hits:
            o.write("\t".join(r) + "\n")

    # 矩阵：每 lead×细胞型；axis=lead 注释基因自身；若 lead 位点不在数据中则 NO_LEAD
    with open(OUT("lead_eQTL_matrix.tsv"), "w", encoding="utf-8") as o:
        o.write("lead_rsid\tlead_gene\tchr\tpos\tin_gene_body\tdataset\tlead_in_data\t"
                "axis_gene_tested\taxis_gene_symbol\taxis_min_p\taxis_beta\taxis_se\t"
                "region_min_p\tregion_min_p_gene\tregion_min_p_gene_symbol\tregion_min_p_beta\n")
        for rsid, lchrom, lpos, gsymbol, gensg in LEADS:
            gb = GENE_BODY.get(gsymbol, ("", 0, 0))
            inbody = "yes" if (gb[0] == lchrom and gb[1] <= lpos <= gb[2]) else \
                     ("near_%dkb" % (min(abs(lpos - gb[1]), abs(lpos - gb[2])) // 1000) if gb[0] == lchrom else "no")
            for dslabel in sorted(set(QTD2LABEL.values())):
                sub = [r for r in all_hits if r[0] == rsid and r[4] == dslabel]
                lead_in = "yes" if sub and sub[0][6] != "NO_LEAD_SNP_IN_REGION" else "no"
                axis = [r for r in sub if r[7] == gsymbol]
                # region min（同 dataset 同染色体）
                reg = [r for r in region_top if r[0] == lchrom and r[1] == dslabel]
                if reg:
                    rp = float(reg[0][6]); rg = reg[0][3]; rgs = sym(reg[0][3]); rb = reg[0][7]
                else:
                    rp = ""; rg = ""; rgs = ""; rb = ""
                if axis:
                    best = min(axis, key=lambda r: float(r[8]))
                    ap, ab, ase = best[8], best[9], best[10]
                else:
                    ap, ab, ase = "", "", ""
                if lead_in == "no":
                    ap = "LEAD_ABSENT"
                o.write("\t".join([rsid, gsymbol, lchrom, str(lpos), inbody, dslabel, lead_in,
                                   "yes" if axis else "no", gsymbol if axis else "",
                                   str(ap), str(ab), str(ase),
                                   str(rp), rg, rgs, str(rb)]) + "\n")

    with open(OUT("lead_region_top.tsv"), "w", encoding="utf-8") as o:
        o.write("chr\tdataset\tdataset_qtd\tgene_id\tgene_symbol\tpos\tpvalue\tbeta\tse\trsid\tvariant\tLD_r2\n")
        for r in region_top:
            o.write("\t".join(r) + "\n")

    print(f"[OUT] {OUT}/lead_hits_all.tsv ({len(all_hits)} rows)")
    print(f"[OUT] {OUT}/lead_eQTL_matrix.tsv")
    print(f"[OUT] {OUT}/lead_region_top.tsv ({len(region_top)} rows)")

if __name__ == "__main__":
    main()
