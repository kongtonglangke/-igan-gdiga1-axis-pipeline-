# -*- coding: utf-8 -*-
"""
step3_power_simulation.py  —  M2b 共定位的**功效分析**（第二轮新增）

背景（审稿意见「共定位阴性结论缺乏功效分析支撑」）
--------------------------------------------------------------
论文原句：25 个 gene×context 共定位中 14 个可算，`none supported a shared
causal variant (maximum PPH4 = 0.031)`，并据此在讨论中断言 axis 与 IgAN
易感性 "genetically separable"。审稿人指出：在弱 eQTL / 弱 GWAS 信息下，
**即便 H4 为真，PPH4 也可能很低**；把「PPH4 低」等同于「无共享因果结构」
是方法学错误。本脚本用**基于模拟的功效分析**量化这一点的严重程度。

方法
----
- 复用论文自己的 `coloc_abf()`（从同目录 step1_coloc_abf.py 动态导入，
  保证与正文结果**完全同一实现**；先验 p1=p2=1e-4, p12=1e-5, W=0.04）。
- 用真实逐位点数据：eQTL Catalogue region_cache（5 个语境）与
  GCST90018866 (IgAN GWAS) 的 beta/se；等位基因方向与前文一致地对齐。
- 模拟 H4 为真：在 eQTL 顶变异处放**单一共享因果变异**，
  eQTL 真效应固定为观测顶效应（即"eQTL 信号是我们的实测信号"），
  GWAS 真效应按网格 δ 取值；两侧各自加抽样噪声 N(0, se_i^2)；
  用 coloc_abf 重算 PPH4。功效 = 重复中 PPH4 ≥ 0.8 的比例。

输出
----
  results/M2b_coloc_power/power_grid.tsv     每个语境的功率网格
  results/M2b_coloc_power/observed_vs_power.tsv  实测 vs 所需效应
  results/M2b_coloc_power/FigureS8_coloc_power.png/.pdf
"""
import os
import sys
import math
import json
import importlib.util

import numpy as np
import pandas as pd

# --- 路径解析：沿用复现包 _paths ---------------------------------------
_C = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _C if os.path.isfile(os.path.join(_C, "_paths.py")) else os.path.dirname(_C))
import _paths  # noqa: E402

RAW = _paths.at("00_rawdata")
IGAN = RAW("GCST90018866", "GCST90018866.h.tsv.gz")
EQTL_CACHE = RAW("eQTL_Catalogue", "region_cache")
OUT = _paths.at("阶段4.5_复现包", "reproducibility", "results", "M2b_coloc_power")
os.makedirs(OUT, exist_ok=True)

# --- 复用 step1 的 coloc_abf（保证同实现）------------------------------
_spec = importlib.util.spec_from_file_location("step1_coloc_abf", os.path.join(_C, "step1_coloc_abf.py"))
step1 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(step1)
coloc_abf = step1.coloc_abf
load_region_cache = step1.load_region_cache
harmonize_beta = step1.harmonize_beta

RNG = np.random.default_rng(20261008)   # 固定种子，可复现

# C1GALT1 在 5 个语境都有 chr7 数据；这是本论文主结论所在
CONTEXTS = [
    ("QTD000608", "OneK1K_B_naive",       "C1GALT1", "ENSG00000106392", "7", 7_229_549),
    ("QTD000607", "OneK1K_B_memory",      "C1GALT1", "ENSG00000106392", "7", 7_229_549),
    ("QTD000606", "OneK1K_B_intermediate", "C1GALT1", "ENSG00000106392", "7", 7_229_549),
    ("QTD000356", "GTEx_v10_blood",       "C1GALT1", "ENSG00000106392", "7", 7_229_549),
    ("QTD000261", "GTEx_v10_kidney_cortex", "C1GALT1", "ENSG00000106392", "7", 7_229_549),
]

# GWAS 真效应网格（每等位基因 log-OR）
DELTAS = [0.0, 0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.20]
N_REP = 4000


def build_igan_index(positions_lo, positions_hi, chrom):
    """只保留 chr7 且落在 [lo,hi] 的行，省内存。"""
    import gzip
    idx = {}
    with gzip.open(IGAN, "rt", encoding="utf-8", errors="replace") as f:
        f.readline()
        for line in f:
            fld = line.rstrip("\n").split("\t")
            if len(fld) < 12:
                continue
            if fld[0] != chrom:
                continue
            try:
                pos = int(fld[1])
            except Exception:
                continue
            if pos < positions_lo or pos > positions_hi:
                continue
            rsid = fld[11]
            if not rsid or rsid == "NA":
                continue
            try:
                # 与 step1.build_igan_index 完全同一元组布局：
                # (chr, pos, beta, se, p, EA, OA) —— harmonize_beta 依赖该顺序
                idx[rsid] = (fld[0], int(fld[1]), float(fld[4]), float(fld[5]),
                             float(fld[7]), fld[2], fld[3])
            except Exception:
                continue
    return idx


def vectorized_coloc(beta1, se1, beta2, se2, p1=1e-4, p2=1e-4, p12=1e-5, W=0.04):
    """与 step1.coloc_abf 同一公式，但按行向量化（每行一个重复）。"""
    var1 = se1 * se1
    var2 = se2 * se2
    r1 = W / (W + var1)
    r2 = W / (W + var2)
    l1 = 0.5 * np.log1p(-r1) + 0.5 * (beta1 ** 2 / var1) * r1
    l2 = 0.5 * np.log1p(-r2) + 0.5 * (beta2 ** 2 / var2) * r2

    def lse(a):
        m = a.max(axis=1, keepdims=True)
        return (m + np.log(np.exp(a - m).sum(axis=1, keepdims=True))).squeeze(1)

    lS1 = lse(l1)
    lS2 = lse(l2)
    lS12 = lse(l1 + l2)
    lS1S2 = lS1 + lS2
    delta = np.clip(lS12 - lS1S2, None, 0.0)
    with np.errstate(divide="ignore"):
        lH3 = np.where(delta > -700, lS1S2 + np.log1p(-np.exp(delta)), lS1S2)
    lpost = np.column_stack([
        np.zeros_like(lS1),
        math.log(p1) + lS1,
        math.log(p2) + lS2,
        math.log(p1) + math.log(p2) + lH3,
        math.log(p12) + lS12,
    ])
    m = lpost.max(axis=1, keepdims=True)
    e = np.exp(lpost - m)
    PPH = e / e.sum(axis=1, keepdims=True)
    return PPH[:, 4]


def main():
    print("==== 共定位功效分析 (M2b power simulation) ====")
    igan = build_igan_index(7_000_000, 7_400_000, "7")
    print("  IgAN chr7 window index:", len(igan))

    grid_rows = []
    obs_rows = []
    for qtd, ctx, sym, ensg, chr_, tss in CONTEXTS:
        df_e = load_region_cache(qtd, chr_, ensg, tss, flank=1_000_000)
        if df_e is None or len(df_e) == 0:
            continue
        recs = [igan.get(r) for r in df_e["rsid"]]
        aligned = [harmonize_beta(a, rec) if rec is not None else (np.nan, "no_igan")
                   for a, rec in zip(df_e["alt"], recs)]
        df_e["b2"] = [x[0] for x in aligned]
        df_e["dir"] = [x[1] for x in aligned]
        df_e["se2"] = [rec[3] if rec is not None else np.nan for rec in recs]
        d = df_e.dropna(subset=["b2", "se2", "beta", "se"])
        d = d[d["dir"].isin(["same", "flip"])]
        b1 = d["beta"].to_numpy(float)
        s1 = d["se"].to_numpy(float)
        b2 = d["b2"].to_numpy(float)
        s2 = d["se2"].to_numpy(float)
        n = len(b1)
        if n < 5:
            continue
        k = int(np.nanargmax(np.abs(b1 / s1)))          # eQTL 顶变异 = 拟定的共享因果变异
        beta1_obs, z1 = b1[k], b1[k] / s1[k]
        p1v = 2 * (1 - 0.5 * (1 + math.erf(abs(z1) / math.sqrt(2))))

        # 实测（不模拟）
        obs = coloc_abf(b1, s1, b2, s2)
        obs_rows.append(dict(context=ctx, n_SNP=n, causal_rsid=d["rsid"].iloc[k],
                             eqtl_beta=beta1_obs, eqtl_se=s1[k], eqtl_z=z1, eqtl_p=p1v,
                             gwas_beta_at_causal=b2[k], gwas_se_at_causal=s2[k],
                             PPH1=obs["PPH1"], PPH3=obs["PPH3"], PPH4=obs["PPH4"]))

        for delta in DELTAS:
            zt1 = beta1_obs / s1[k]
            zt2 = delta / s2[k]
            R = N_REP
            # 模拟：只有 k 位点有真效应
            Z1 = RNG.standard_normal((R, n))
            Z2 = RNG.standard_normal((R, n))
            Z1[:, k] += zt1
            Z2[:, k] += zt2
            B1 = Z1 * s1
            B2 = Z2 * s2
            pph4 = vectorized_coloc(B1, s1[None, :], B2, s2[None, :])
            grid_rows.append(dict(context=ctx, n_SNP=n, causal_rsid=d["rsid"].iloc[k],
                                  eqtl_z=zt1, gwas_delta=delta, gwas_OR=math.exp(delta),
                                  gwas_z_at_causal=zt2,
                                  power_pph4_ge_0_8=float((pph4 >= 0.8).mean()),
                                  power_pph4_ge_0_5=float((pph4 >= 0.5).mean()),
                                  median_pph4=float(np.median(pph4))))
            print(f"  {ctx:24s} n={n:5d} delta={delta:5.3f} "
                  f"pow(0.8)={(pph4>=0.8).mean():6.3f} pow(0.5)={(pph4>=0.5).mean():6.3f}")

    pd.DataFrame(grid_rows).to_csv(os.path.join(OUT, "power_grid.tsv"), sep="\t", index=False)
    pd.DataFrame(obs_rows).to_csv(os.path.join(OUT, "observed_vs_power.tsv"), sep="\t", index=False)
    print("saved power_grid.tsv / observed_vs_power.tsv")

    # ---- 图：功率曲线 ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        g = pd.DataFrame(grid_rows)
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), dpi=300)
        order = ["OneK1K_B_naive", "OneK1K_B_memory", "OneK1K_B_intermediate",
                 "GTEx_v10_blood", "GTEx_v10_kidney_cortex"]
        for ctx in order:
            sub = g[g["context"] == ctx].sort_values("gwas_delta")
            if sub.empty:
                continue
            axes[0].plot(sub["gwas_OR"], sub["power_pph4_ge_0_8"], marker="o", label=ctx)
            axes[1].plot(sub["gwas_OR"], sub["power_pph4_ge_0_5"], marker="o", label=ctx)
        for ax, ttl in zip(axes, ["Power at PPH4 >= 0.8", "Power at PPH4 >= 0.5"]):
            ax.set_xlabel("True shared-variant effect (OR per allele)")
            ax.set_ylabel("Power")
            ax.set_ylim(-0.03, 1.03)
            ax.grid(alpha=.3)
            ax.set_title(ttl)
            ax.axhline(0.8, ls=":", c="grey", lw=1)
        axes[0].legend(fontsize=7, loc="upper left")
        fig.suptitle("Colocalisation power at the C1GALT1 locus (simulation under shared-variant H4)")
        fig.tight_layout()
        for ext in ("png", "pdf", "tiff"):
            fig.savefig(os.path.join(OUT, f"FigureS8_coloc_power.{ext}"))
        print("saved FigureS8_coloc_power.png/.pdf/.tiff")
    except Exception as exc:  # pragma: no cover
        print("plot skipped:", exc)


if __name__ == "__main__":
    main()
