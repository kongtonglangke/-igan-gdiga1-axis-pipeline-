#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
step3_motif_cpg.py — 甲基化敏感调控逻辑：轴位点 CpG × TF motif 交叉
================================================================================
问题：L3 已把 C1GALT1 的 5 个显著 CpG 定性为"转录耦合的基因体甲基化标记"，并把
      C1GALT1C1 启动子 CpG island 高甲基化列为 L6 干预窗的锚点。但"甲基化如何
      影响调控"这一步缺一个**序列层面**的可检验假说：
      哪些 TF 的结合位点与这些 CpG 重叠（即其结合可能被 CpG 甲基化直接阻断）？

流程
----
  (1) 取序列：
      · C1GALT1 5 个显著基因体 CpG（hg19 坐标）± 60 bp
      · C1GALT1C1 启动子 CpG island（hg38 chrX:120629882-120630375，跨 TSS）
  (2) JASPAR 2024 CORE 脊椎动物非冗余 PFM（879 motif）→ log-odds PWM 扫描，
      阈值 = motif 最大可能得分的 80%。
  (3) 标注命中窗口是否**覆盖 CpG 二核苷酸** → 甲基化敏感候选。
  (4) 与 step2 推断出的轴关联调控子取交集，产出 Table S9。

输出：results/sc_regulatory/motif_hits_axis_cpgs.tsv（全量命中）
      results/sc_regulatory/TableS9_methylation_sensitive_TFs.tsv（汇总）
"""

# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import os
import io
import re
import json
import subprocess

import numpy as np
import pandas as pd

BASE = _paths.LEGACY
DATA_DIR = BASE("阶段4.5_复现包", "reproducibility", "data", "sc_regulatory")
RES_DIR = BASE("阶段4.5_复现包", "reproducibility", "results", "sc_regulatory")
os.makedirs(RES_DIR, exist_ok=True)

JASPAR_URL = ("https://jaspar.genereg.net/download/data/2024/CORE/"
              "JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt")
JASPAR_LOCAL = os.path.join(DATA_DIR, "JASPAR2024_CORE_vertebrates_nr.txt")
UCSC = "https://api.genome.ucsc.edu/getData/sequence?genome={g};chrom={c};start={s};end={e}"
FLANK = 60
REL_THRESH = 0.80

# C1GALT1 五个显著基因体 CpG（GoDMC/450K，hg19/GRCh37 坐标）
CPGS = [
    ("cg19603390", "chr7", 7222222),
    ("cg19473623", "chr7", 7224869),
    ("cg17994788", "chr7", 7261594),
    ("cg04827551", "chr7", 7268805),
    ("cg16101574", "chr7", 7291514),
]
# C1GALT1C1 启动子 CpG island（hg38）
C1GALT1C1_ISLAND = ("C1GALT1C1_promoter_CpG_island", "chrX", 120629882, 120630375, "hg38")


def fetch(url):
    r = subprocess.run(["curl", "-sSL", "--ssl-no-revoke", "-m", "60", url],
                       capture_output=True, text=True)
    return r.stdout


def get_seq(genome, chrom, start, end):
    out = fetch(UCSC.format(g=genome, c=chrom, s=start, e=end))
    return json.loads(out).get("dna", "").upper()


def load_jaspar():
    """返回 {motif_id: (tf_name, pfm)}。JASPAR jaspar 格式表头为 '>MA0006.2\tArnt'。"""
    if not os.path.exists(JASPAR_LOCAL):
        txt = fetch(JASPAR_URL)
        if not txt.strip().startswith(">"):
            raise RuntimeError("JASPAR 下载失败")
        open(JASPAR_LOCAL, "w", encoding="utf-8").write(txt)
    motifs = {}
    mid, tfname, rows = None, None, {}

    def flush():
        if mid and len(rows) == 4:
            motifs[mid] = (tfname, np.array([rows[b] for b in "ACGT"], dtype=float))

    for line in open(JASPAR_LOCAL, encoding="utf-8"):
        line = line.rstrip("\n")
        if line.startswith(">"):
            flush()
            parts = line[1:].split("\t")
            mid = parts[0].strip()
            tfname = parts[1].strip() if len(parts) > 1 and parts[1].strip() else mid
            rows = {}
        elif line and line[0] in "ACGT":
            nums = re.findall(r"-?\d+(?:\.\d+)?", line[1:])
            rows[line[0]] = [float(x) for x in nums]
    flush()
    return motifs


def pwm(pfm):
    n = pfm.sum(axis=0, keepdims=True)
    p = (pfm + 0.25 * (n + 1)) / (n + 1)
    return np.log2(p / 0.25)


def scan(seq, motifs, rel_thresh=REL_THRESH):
    idx = {b: i for i, b in enumerate("ACGT")}
    codes = np.array([idx.get(c, -1) for c in seq])
    hits = []
    for mid, (tfname, pfm) in motifs.items():
        W = pwm(pfm)
        L = W.shape[1]
        if L > len(seq):
            continue
        max_score = W.max(axis=0).sum()
        thresh = rel_thresh * max_score
        # 滑动窗口得分
        for s in range(len(seq) - L + 1):
            win = codes[s:s + L]
            if (win < 0).any():
                continue
            sc = W[win, np.arange(L)].sum()
            if sc >= thresh:
                hits.append({"motif": mid, "tf": tfname,
                             "start": s, "end": s + L, "score": sc,
                             "rel": sc / max_score})
    return hits


def cpg_positions(seq):
    return [i for i in range(len(seq) - 1) if seq[i] == "C" and seq[i + 1] == "G"]


def main():
    print("=" * 74)
    print("step3 甲基化敏感调控逻辑：轴位点 CpG × JASPAR motif")
    print("=" * 74)
    motifs = load_jaspar()
    print(f"JASPAR motif 数 = {len(motifs)}")

    regions = []
    for cpg, chrom, pos in CPGS:
        seq = get_seq("hg19", chrom, pos - FLANK, pos + FLANK + 1)
        regions.append({"region": cpg, "group": "C1GALT1_gene_body_CpG",
                        "genome": "hg19", "chrom": chrom, "start": pos - FLANK,
                        "seq": seq, "cpg_offset": FLANK})
        print(f"  {cpg:12s} {chrom}:{pos}  长度 {len(seq)}  CpG@ {FLANK}")
    nm, c, s, e, g = C1GALT1C1_ISLAND
    seq = get_seq(g, c, s, e)
    regions.append({"region": nm, "group": "C1GALT1C1_promoter_island",
                    "genome": g, "chrom": c, "start": s, "seq": seq, "cpg_offset": None})
    print(f"  {nm}  {c}:{s}-{e}  长度 {len(seq)}  含 {len(cpg_positions(seq))} 个 CpG")

    all_hits = []
    for r in regions:
        hs = scan(r["seq"], motifs)
        cps = set(cpg_positions(r["seq"]))
        for h in hs:
            covered = [p for p in cps if h["start"] <= p <= h["end"] - 1]
            all_hits.append({**{k: r[k] for k in ("region", "group", "genome", "chrom", "start")},
                             "motif": h["motif"], "tf": h["tf"], "rel_score": round(h["rel"], 4),
                             "motif_start_in_window": h["start"], "motif_end_in_window": h["end"],
                             "covers_CpG": bool(covered),
                             "n_CpG_in_motif": len(covered)})
        print(f"  [{r['region']}] motif 命中 {len(hs)}，覆盖 CpG 的 {sum(1 for x in all_hits if x['region']==r['region'] and x['covers_CpG'])}")

    df = pd.DataFrame(all_hits)
    df.to_csv(os.path.join(RES_DIR, "motif_hits_axis_cpgs.tsv"), sep="\t", index=False)

    # ---- 汇总：甲基化敏感候选（覆盖 CpG 的 motif，按区域/去重 TF） ----
    ms = df[df["covers_CpG"]].copy()
    summ = (ms.groupby(["region", "group", "tf"])
              .agg(n_hits=("motif", "size"), best_rel=("rel_score", "max"),
                   max_CpG=("n_CpG_in_motif", "max"))
              .reset_index().sort_values(["group", "best_rel"], ascending=[True, False]))
    summ.to_csv(os.path.join(RES_DIR, "TableS9_methylation_sensitive_TFs.tsv"),
                sep="\t", index=False)
    print(f"\n甲基化敏感候选总数（区域×TF）= {len(summ)}，涉及 TF {summ['tf'].nunique()} 个")

    # ---- 与 step2 调控子交集 ----
    # 主判据：q_emp < 0.10 —— step2 的「按自身靶数匹配零分布」经验 P 经 BH-FDR 校正
    #         （唯一被本文认可的显著性口径）
    # 参考判据：emp_p < 0.05（未校正名义）或 top40（宽松参考）
    comp_path = os.path.join(RES_DIR, "regulon_axis_composite.tsv")
    inter, inter_broad = [], []
    if os.path.exists(comp_path):
        comp = pd.read_csv(comp_path, sep="\t")
        hit_tfs = set(summ["tf"])
        sig = comp[comp["q_emp"] < 0.10]
        inter = sorted((set(sig["tf"]) & hit_tfs))
        broad = set(comp[comp["emp_p"] < 0.05]["tf"]) | set(comp.head(40)["tf"])
        inter_broad = sorted(broad & hit_tfs)
        print(f"\n[交集·q_emp<0.10] n={len(inter)}: {inter}")
        print(f"[交集·emp_p<0.05|top40] n={len(inter_broad)}: {inter_broad[:25]}")
        n_emp05 = int((comp["emp_p"] < 0.05).sum())
        n_q10 = int((comp["q_emp"] < 0.10).sum())
    else:
        n_emp05 = n_q10 = -1
        print("\n[警告] 未找到 regulon_axis_composite.tsv —— 请先运行 step2 再运行 step3")

    with open(os.path.join(RES_DIR, "step3_summary.json"), "w", encoding="utf-8") as f:
        json.dump({"n_motifs": len(motifs), "n_hits": int(len(df)),
                   "n_methylation_sensitive_raw_hits": int(len(ms)),
                   "n_region_tf_rows": int(len(summ)),
                   "n_methylation_sensitive_hits": int(len(summ)),
                   "n_candidate_TFs": int(summ["tf"].nunique()),
                   "step2_n_passed_emp05": n_emp05,
                   "step2_n_passed_qemp010": n_q10,
                   "intersect_with_step2": inter,
                   "intersect_with_step2_broad": inter_broad[:60]}, f, indent=2, ensure_ascii=False)
    print("\n[done] step3")


if __name__ == "__main__":
    main()
