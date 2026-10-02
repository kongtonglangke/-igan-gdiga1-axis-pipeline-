# -*- coding: utf-8 -*-
"""阶段2 主证据1 前奏：4 lead Gd-IgA1 位点 × OneK1K B 细胞型 + GTEx v10 blood 的 cis-eQTL 查询

方法：用 pysam.TabixFile 通过 https 远程打开 eQTL Catalogue r8 .all.tsv.gz.tbi（permuted 只给 top；
all.tsv.gz 含全 nominal 关联，tabix 索引存在），按 lead 位点 ±500kb 区域 fetch，输出 top cis-eQTL。
lead 位点（GRCh38 坐标 from Kiryluk/Wang 文献核验表）：
  rs13226913  chr7:7207215  C1GALT1
  rs5910940   X:120624955  C1GALT1C1
  rs10238682  chr7:7215386  C1GALT1
  rs7856182   chr9:98871030 GALNT12
注：rs5910940 X 染色体 — eQTL Cat r8 OneK1K 不收 X 染色体 eQTL；C1GALT1C1 在 B 细胞 no-hit 已知，符合预期。

依赖：pysam、requests（已装 pysam venv 即可）
输出：阶段2/主证据1_lead_eQTL/lead_eQTL_top.tsv + 详细 log
"""
import os, gzip, sys, time
import pysam

LEADS = [
    # (rsID, chrom, pos_hg38, gene)
    ("rs13226913", "7",  7207215,  "C1GALT1"),
    ("rs5910940",  "X",  120624955,"C1GALT1C1"),
    ("rs10238682", "7",  7215386,  "C1GALT1"),
    ("rs7856182",  "9",  98871030, "GALNT12"),
]

# (label, qts, qtd)
DATASETS = [
    ("OneK1K_B_naive",         "QTS000038", "QTD000608"),
    ("OneK1K_B_memory",        "QTS000038", "QTD000607"),
    ("OneK1K_B_intermediate",  "QTS000038", "QTD000606"),
    ("GTEx_v10_blood",         "QTS000015", "QTD000356"),
    ("GTEx_v10_kidney_cortex", "QTS000015", "QTD000261"),
]

BASE = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats"

def open_tabix(label, qts, qtd):
    url = f"{BASE}/{qts}/{qtd}/{qtd}.all.tsv.gz"
    print(f"  open {url}", flush=True)
    return pysam.TabixFile(url)

def query(tabix, chrom, pos, flank=500_000):
    region = f"{chrom}:{max(1, pos - flank)}-{pos + flank}"
    try:
        rows = list(tabix.fetch(region, parser=pysam.asTuple()))
    except Exception as e:
        print(f"  WARN fetch {region}: {e}")
        return []
    return rows

# 1) 列结构先抓一个 lead × 一个 dataset 看列
print("[STEP1] sample fetch — rs13226913 × OneK1K B_naive")
t0 = time.time()
tb = open_tabix(*DATASETS[0])
sample = query(tb, "7", 7207215)
print(f"  {len(sample)} rows in ±500kb of rs13226913 (B_naive); first row:", sample[0] if sample else "EMPTY")
print(f"  header by column call: {tb.header}")
tb.close()
print(f"  sample roundtrip {time.time()-t0:.1f}s")

# 2) 完整查询 4 lead × 5 datasets
print("\n[STEP2] full 4 lead × 5 dataset query")
out = []
log = [f"opened at {time.strftime('%Y-%m-%d %H:%M:%S')}, {len(LEADS)} leads × {len(DATASETS)} datasets"]
for ds_label, qts, qtd in DATASETS:
    tb = open_tabix(ds_label, qts, qtd)
    for rsid, chrom, pos, gene in LEADS:
        rows = query(tb, chrom, pos)
        # 关联表列（来自 eQTL Cat Columns.md）：0 molecular_trait_id 1 chromosome 2 position
        # 3 ref 4 alt 5 variant 6 rsids 7 tss_distance 8 ma_samples 9 pvalue 10 an 11 beta
        # 12 se 13 p_perm 14 p_beta 15 q_beta 16 q_perm 17 q_value
        keep = []
        for r in rows:
            try:
                if int(r[2]) == pos and rsid in r[6].split(","):
                    keep.append(r)
            except Exception:
                continue
        if not keep:
            out.append([rsid, gene, ds_label, "", "", "no_lead_eQTL", "", "", ""])
            continue
        # 选 pvalue 最小
        def pval(r):
            try: return float(r[9])
            except Exception: return 1.0
        best = min(keep, key=pval)
        out.append([rsid, gene, ds_label, best[0], best[1] + ":" + best[2] + ":" + best[3] + ":" + best[4],
                    best[9], best[11], best[12]])
    tb.close()

# 3) 输出
os.makedirs("阶段2/主证据1_lead_eQTL", exist_ok=True)
out_path = "阶段2/主证据1_lead_eQTL/lead_eQTL_top.tsv"
with open(out_path, "w", encoding="utf-8") as f:
    f.write("rsid\tgene\tdataset\tmolecular_trait_id\ttop_variant\tpvalue_nominal\tbeta\tse\n")
    for r in out:
        f.write("\t".join(r) + "\n")
print(f"\n[OUT] {out_path}")
print("\n=== 概要 ===")
for r in out:
    print("  ", r)
