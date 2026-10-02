#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
step1_bcell_eqtl_heterogeneity.py — B2：C1GALT1 lead 变异的 B 亚型 eQTL 异质性检验
================================================================================
问题：现稿只说 C1GALT1 lead 变异在 B 细胞有 cis-eQTL，但**未检验效应是否随 B 细胞
      分化状态改变**。若效应只存在于 naive B（而记忆/中间态没有），则说明该遗传效应
      被 B 细胞状态"门控"，为"上游调控层 L7"提供遗传学呼应。

数据：eQTL Catalogue r8 OneK1K（Yazar et al. Science 2022）B 亚型 nominal sumstats
      缓存于 00_rawdata/eQTL_Catalogue/region_cache/QTD0006{08,07,06}.7.tsv
        QTD000608 = B_naive        （OneK1K 14 细胞型之一）
        QTD000607 = B_memory
        QTD000606 = B_intermediate （unswitched / marginal-zone-like）
      列（1-based）：1 molecular_trait_id(ENSG) 2 chromosome 3 position 4 ref 5 alt
                     6 variant 7 ma_samples 8 af 9 pvalue 10 beta 11 se
                     12 variant_type 13 n 14 an 15 p_perm 16/17 trait 18 NA 19 rsids
      （注：缓存文件首行为下载截断的残行，不含 ENSG，按 $1==ENSG 过滤即可跳过）

统计：
  · 固定效应合并 beta；Cochran's Q 异质性（df=k-1=2）；I²
  · 成对 z 检验（naive vs memory / naive vs intermediate / memory vs intermediate）
  · 加权元回归对比：naive vs 抗原经历态（memory+intermediate）
  · 依发育序的加权趋势检验（naive=0, intermediate=1, memory=2）

输出：results/B2_bcell_eqtl/TableS10_Bcell_eQTL_heterogeneity.tsv
      results/B2_bcell_eqtl/bcell_eqtl_effects.tsv
      results/B2_bcell_eqtl/step1_summary.json
      results/B2_bcell_eqtl/Figure_S7_Bcell_eQTL_heterogeneity.{png,pdf}
"""

# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import os
import json

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = _paths.LEGACY
CACHE = BASE("00_rawdata", "eQTL_Catalogue", "region_cache")
RES = BASE("阶段4.5_复现包", "reproducibility", "results", "B2_bcell_eqtl")
os.makedirs(RES, exist_ok=True)

# (rsid, chrom_cache_tag, pos_hg38, trait_ensg, gene_symbol)
LEADS = [
    ("rs13226913", "7", 7207215, "ENSG00000106392", "C1GALT1"),
    ("rs10238682", "7", 7215386, "ENSG00000106392", "C1GALT1"),
]
# OneK1K B 亚型（发育序：naive → intermediate → memory）
SUBSETS = [
    ("B_naive", "QTD000608", 0),
    ("B_intermediate", "QTD000606", 1),
    ("B_memory", "QTD000607", 2),
]
STAGE_ORDER = {"B_naive": 0, "B_intermediate": 1, "B_memory": 2}


def read_effect(qtd, chrom_tag, pos, ensg, rsid):
    """从 region_cache 抽 (beta, se, p, af, n) —— 只在匹配行上解析，规避残行。"""
    path = os.path.join(CACHE, f"{qtd}.{chrom_tag}.tsv")
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 19:
                continue
            if p[0] == ensg and p[2] == str(pos):
                try:
                    return {"beta": float(p[9]), "se": float(p[10]), "p": float(p[8]),
                            "af": float(p[7]), "n": int(float(p[6])), "rsid": p[18]}
                except ValueError:
                    return None
    return None


def fe_meta(b, se):
    """固定效应合并 + Cochran's Q + I²。"""
    w = 1.0 / np.asarray(se) ** 2
    b = np.asarray(b)
    b_fe = float((w * b).sum() / w.sum())
    se_fe = float(np.sqrt(1.0 / w.sum()))
    Q = float((w * (b - b_fe) ** 2).sum())
    dfree = len(b) - 1
    p_Q = float(stats.chi2.sf(Q, dfree)) if dfree > 0 else np.nan
    I2 = float(max(0.0, (Q - dfree) / Q) * 100) if Q > 0 else 0.0
    return b_fe, se_fe, Q, dfree, p_Q, I2, w


def wls_slope(x, b, se):
    """加权最小二乘 b = a + c*x；返回 c 及其 SE / p（固定效应假定）。"""
    X = np.column_stack([np.ones(len(x)), np.asarray(x, float)])
    w = 1.0 / np.asarray(se, float) ** 2
    W = np.diag(w)
    XtWX = X.T @ W @ X
    cov = np.linalg.inv(XtWX)
    coef = cov @ (X.T @ W @ np.asarray(b, float))
    c = float(coef[1])
    se_c = float(np.sqrt(cov[1, 1]))
    z = c / se_c
    return c, se_c, float(z), float(2 * stats.norm.sf(abs(z)))


def main():
    print("=" * 74)
    print("B2：C1GALT1 lead × OneK1K B 亚型 eQTL 异质性")
    print("=" * 74)

    rows, tests = [], []
    for rsid, chrom_tag, pos, ensg, sym in LEADS:
        est = []
        for subset, qtd, _ in SUBSETS:
            r = read_effect(qtd, chrom_tag, pos, ensg, rsid)
            if r is None:
                print(f"  [!] {rsid} × {subset}: 无记录")
                continue
            est.append({"lead_rsid": rsid, "gene": sym, "subset": subset, **r})
            rows.append({"lead_rsid": rsid, "gene": sym, "subset": subset,
                         "beta": r["beta"], "se": r["se"], "p": r["p"],
                         "af": r["af"], "ma_samples": r["n"]})
            print(f"  {rsid} {sym} {subset:15s} beta={r['beta']:+.4f} se={r['se']:.4f} p={r['p']:.3g}")

        if len(est) < 2:
            continue
        df = pd.DataFrame(est)
        # 按发育序排（naive → intermediate → memory）
        df = df.sort_values("subset", key=lambda c: c.map(STAGE_ORDER))
        b, se = df["beta"].values, df["se"].values
        b_fe, se_fe, Q, dfree, p_Q, I2, w = fe_meta(b, se)

        def rec(name, val, se_v=None, p=None, note=""):
            tests.append({"lead_rsid": rsid, "gene": sym, "test": name,
                          "stat": val, "se": se_v, "p": p, "note": note})

        rec("fixed_effect_beta", b_fe, se_fe,
            float(2 * stats.norm.sf(abs(b_fe / se_fe))), "inverse-variance pooling of 3 subsets")
        rec("cochran_Q", Q, None, p_Q, f"df={dfree}; I2={I2:.1f}%")

        # 成对 z 检验
        idx = {s: i for i, s in enumerate(df["subset"])}
        for a, bb in [("B_naive", "B_intermediate"), ("B_naive", "B_memory"),
                      ("B_intermediate", "B_memory")]:
            if a in idx and bb in idx:
                d = b[idx[a]] - b[idx[bb]]
                sd = float(np.sqrt(se[idx[a]] ** 2 + se[idx[bb]] ** 2))
                z = d / sd
                rec(f"pairwise_{a}_vs_{bb}", float(d), sd,
                    float(2 * stats.norm.sf(abs(z))), "difference in beta")

        # 元回归：naive vs 抗原经历态
        x = np.array([0.0 if s == "B_naive" else 1.0 for s in df["subset"]])
        c, se_c, z, p = wls_slope(x, b, se)
        rec("metareg_naive_vs_experienced", c, se_c, p,
            "x=0 naive, x=1 memory|intermediate")

        # 趋势检验（发育序）
        xr = np.array([STAGE_ORDER[s] for s in df["subset"]], float)
        c, se_c, z, p = wls_slope(xr, b, se)
        rec("trend_per_stage", c, se_c, p, "rank naive<intermediate<memory (WLS slope)")

    eff = pd.DataFrame(rows)
    eff.to_csv(os.path.join(RES, "bcell_eqtl_effects.tsv"), sep="\t", index=False)
    td = pd.DataFrame(tests)
    td.to_csv(os.path.join(RES, "TableS10_Bcell_eQTL_heterogeneity.tsv"), sep="\t", index=False)

    print("\n--- 检验结果 ---")
    for _, t in td.iterrows():
        pv = "" if pd.isna(t["p"]) else f" p={t['p']:.3g}"
        print(f"  {t['lead_rsid']:11s} {t['test']:34s} stat={t['stat']:+.4f}{pv}")

    # ---- 森林图 ----
    fig, axes = plt.subplots(1, len(LEADS), figsize=(9.4, 3.2), sharey=True)
    if len(LEADS) == 1:
        axes = [axes]
    for ax, (rsid, chrom_tag, pos, ensg, sym) in zip(axes, LEADS):
        d = eff[eff["lead_rsid"] == rsid].copy()
        d = d.sort_values("subset", key=lambda c: c.map(STAGE_ORDER))
        y = np.arange(len(d))[::-1]
        ax.errorbar(d["beta"], y, xerr=1.96 * d["se"], fmt="o", color="#b2182b",
                    ecolor="#666666", capsize=3, ms=6)
        ax.axvline(0, color="grey", lw=0.8, ls="--")
        ax.set_yticks(y)
        ax.set_yticklabels([s.replace("B_", "") for s in d["subset"]], fontsize=9)
        ax.set_xlabel(f"cis-eQTL β on {sym} (per allele)", fontsize=9)
        for yi, (bb, sp) in zip(y, zip(d["beta"], d["p"])):
            ax.text(bb, yi + 0.22, f"p={sp:.1g}", fontsize=7, ha="center", color="#333333")
        ax.set_title(rsid, fontsize=10)
        ax.tick_params(labelsize=8)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    fig.suptitle("C1GALT1 cis-eQTL effect by B-cell state (OneK1K)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(RES, f"Figure_S7_Bcell_eQTL_heterogeneity.{ext}"), dpi=300)
    plt.close(fig)

    summary = {
        "n_leads": len(LEADS),
        "subsets": [s for s, _, _ in SUBSETS],
        "effects": rows,
        "tests": tests,
    }
    with open(os.path.join(RES, "step1_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print("\n[done] B2")


if __name__ == "__main__":
    main()
