# -*- coding: utf-8 -*-
"""
阶段1 M2a：多细胞语境 cis-eQTL 提取与比较（eQTL Catalogue r8 通道）
数据源：eQTL Catalogue r8 permuted（每基因×每数据集 最显著 cis-eQTL）
  OneK1K  B_naive(B)/B_memory(B)/B_intermediate(B) + GTEx_v10 全血/肾皮质
  （注：OneK1K plasma 未收录于 eQTL Cat r8；eQTLGen 不在 eQTL Cat → 另源补）
输入：00_rawdata/eQTL_Catalogue/*.permuted.tsv.gz
输出：<REPO>/data/M2a_eqtl/ 下的 eQTL_5genes_x_context.tsv、
      eQTL_5genes_summary_heat.tsv 与 format_check.log
用法：python screen_eqtlcat_5genes.py
"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import gzip, os, glob, csv

RAW = _paths.at("00_rawdata/eQTL_Catalogue")     # 外部原始数据（不随仓库分发，需自行下载）
OUT_DIR = _paths.at("阶段1/M2a_eQTL语境")        # → <REPO>/data/M2a_eqtl
GENES = {
    "ENSG00000106392": "C1GALT1",
    "ENSG00000171155": "C1GALT1C1",
    "ENSG00000143641": "GALNT2",
    "ENSG00000119514": "GALNT12",
    "ENSG00000070731": "ST6GALNAC2",
}
# context 标签（取自 dataset_metadata_r8：QTS/QTD）
CONTEXTS = {
    "OneK1K_B_naive":        ("QTS000038", "QTD000608", "B naive", 876),
    "OneK1K_B_memory":       ("QTS000038", "QTD000607", "B memory", 726),
    "OneK1K_B_intermediate": ("QTS000038", "QTD000606", "B intermediate", 618),
    "GTEx_v10_blood":        ("QTS000015", "QTD000356", "GTEx v10 blood (Whole_Blood)", 853),
    "GTEx_v10_kidney_cortex":("QTS000015", "QTD000261", "GTEx v10 kidney cortex", 106),
}

found_gene_in = {g: [] for g in GENES}     # 记录 ENSG 出现在哪些 context
out_rows = []
log = []
for label, (qts, qtd, desc, n) in CONTEXTS.items():
    fp = RAW(f"{label}.permuted.tsv.gz")
    if not os.path.exists(fp):
        log.append(f"[WARN] {label}: 文件缺失 {fp}")
        continue
    with gzip.open(fp, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        col = {c: i for i, c in enumerate(header)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            ens = p[col["molecular_trait_id"]] if "molecular_trait_id" in col else ""
            # 兼容旧列名
            if not ens.startswith("ENSG") and "gene_id" in col:
                ens = p[col["gene_id"]]
            if ens not in GENES:
                continue
            found_gene_in[ens].append(label)
            rec = {
                "gene": GENES[ens], "ensg": ens,
                "dataset": label, "dataset_desc": desc, "n_individuals": n,
                "top_variant": p[col["variant"]] if "variant" in col else "",
                "chr": p[col["chromosome"]] if "chromosome" in col else "",
                "pos": p[col["position"]] if "position" in col else "",
                "beta": p[col["beta"]] if "beta" in col else "",
                "p_nominal_top": p[col["pvalue"]] if "pvalue" in col else "",
                "p_perm": p[col["p_perm"]] if "p_perm" in col else "",
                "p_beta": p[col["p_beta"]] if "p_beta" in col else "",
            }
            out_rows.append(rec)
    log.append(f"[OK] {label} ({desc}, n={n}): 表头 {len(header)} 列，含分子ID列")
    # 行数统计
    nrow = sum(1 for _ in gzip.open(fp, "rt"))
    log.append(f"     总基因×行数={nrow}")

# ---- 输出 1：5基因×语境 全命中表 ----
os.makedirs(OUT_DIR, exist_ok=True)
out1 = OUT_DIR("eQTL_5genes_x_context.tsv")
with open(out1, "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["gene", "ensg", "dataset", "dataset_desc", "n_individuals",
                "top_variant", "chr", "pos", "beta", "p_nominal_top", "p_perm", "p_beta"])
    w.writerows([[r["gene"], r["ensg"], r["dataset"], r["dataset_desc"], r["n_individuals"],
                  r["top_variant"], r["chr"], r["pos"], r["beta"], r["p_nominal_top"],
                  r["p_perm"], r["p_beta"]] for r in out_rows])

# ---- 输出 2：汇总热度（每 基因×context 是否可测 + 最显著）----
def pf(x):
    try:
        return float(x)
    except Exception:
        return None

summary = []
for ens, sym in GENES.items():
    for label, (qts, qtd, desc, n) in CONTEXTS.items():
        rows = [r for r in out_rows if r["ensg"] == ens and r["dataset"] == label]
        if not rows:
            summary.append([sym, ens, label, desc, "", "", "NA", "NA", "NA", "no-hit"])
            continue
        # 取 p_perm 最小者（若有多个转录本/行）
        best = min(rows, key=lambda r: pf(r["p_perm"]) if pf(r["p_perm"]) is not None else 1)
        sig = "sig" if (pf(best["p_perm"]) is not None and pf(best["p_perm"]) < 0.05) else ("nominal" if pf(best["p_nominal_top"]) is not None and pf(best["p_nominal_top"]) < 0.05 else "ns")
        summary.append([sym, ens, label, desc, best["top_variant"], best["chr"] + ":" + best["pos"],
                        best["beta"], best["p_perm"], best["p_beta"], sig])
out2 = OUT_DIR("eQTL_5genes_summary_heat.tsv")
with open(out2, "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["gene", "ensg", "dataset", "dataset_desc", "top_variant", "chr_pos",
                "beta", "p_perm", "p_beta", "call"])
    w.writerows(summary)

# ---- 输出 3：格式校验日志 ----
for ens, sym in GENES.items():
    log.append(f"[GENE] {sym} {ens}: 出现在 {len(found_gene_in[ens])} 个数据集 -> {found_gene_in[ens]}")
with open(OUT_DIR("format_check.log"), "w", encoding="utf-8") as f:
    f.write("\n".join(log) + "\n")
print("\n".join(log))
print(f"\n[OUT] {out1}\n[OUT] {out2}")
