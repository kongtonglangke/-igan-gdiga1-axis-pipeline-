#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
step2_regulon.py — 上游调控子（regulon）推断及其与轴基因的关联（GSE285335）
================================================================================
问题：谁在 IgAN B 细胞中与 O-糖基化轴基因的表达协同？CollecTRI 未收录轴基因
      为靶点，因此采用**数据驱动的调控子活性—轴表达关联**：
  (1) 用 CollecTRI（OmniPath）TF→靶基因网络（限定在表达面板内，≥5 靶）构建调控子。
      **不设靶数上限**（2026-10-01 起；完整网络下若沿用旧上限会人为排除大调控子）。
  (2) 面板基因按 log1p CPM 跨细胞 z 标准化；TF 活性 = 靶基因 z 均值（ULM 型，
      decoupleR 同源做法）。
  (3) **控制细胞类型混杂**：仅在 B 区室内，按 (供者 × B 亚群) 伪批量 → 跨供者
      Spearman 相关（TF 活性 vs 轴基因表达），n = 供者数。
  (4) **控制靶数混杂（本版核心变更）**：经验零分布由"固定按中位靶数的随机基因集"
      改为**按每个调控子自身靶数匹配**的 1,000 个随机基因集。旧口径下，靶数 8 的
      KLF2 与靶数 429 的 SP1 都被拿去和同一批 15 靶随机集比较，量纲不对齐；
      完整网络的靶数跨度更大（5–895），该问题被放大。
  (5) 显著性：逐 TF 经验 P（emp_p）→ **同一 B 亚群内 BH-FDR**（q_emp）。

输出：results/sc_regulatory/regulon_*.tsv + Figure_S6

依赖：`code/lib/collectri.py` —— CollecTRI 属"不随仓库重分发"的第三方资源，
      缺失时由该模块自动从 OmniPath 获取（并打印与论文所用网络的版本差异告警）。
"""

# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
_sys.path.insert(0, _paths.LIB)                  # 随包分发的辅助模块
from collectri import load as load_collectri  # noqa: E402  （CollecTRI 缺失时自动获取）
# --------------------------------------------------------------------------
import os
import json
import argparse

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import spearmanr, rankdata
from statsmodels.stats.multitest import multipletests

BASE = _paths.LEGACY
DATA_DIR = BASE("阶段4.5_复现包", "reproducibility", "data", "sc_regulatory")
RES_DIR = BASE("阶段4.5_复现包", "reproducibility", "results", "sc_regulatory")
FIG_DIR = _paths.results("sc_regulatory")  # 原目的地为投稿材料 additional_file_1
for d in (RES_DIR, FIG_DIR):
    os.makedirs(d, exist_ok=True)

AXIS = ["C1GALT1", "C1GALT1C1", "GALNT2", "GALNT12", "ST6GALNAC2"]
SUBSETS = ["B_naive", "B_memory", "B_plasma"]
# 展示用：与"细胞因子 → O-糖基化轴 / 表观机制"叙事直接相关的 TF（预先指定，勿随结果调整）
CYTOKINE_TFS = ["STAT1", "STAT3", "STAT6", "NFKB1", "RELA", "NFKB2",
                "IRF1", "IRF4", "IRF8", "JUN", "JUNB", "FOS",
                "PRDM1", "PAX5", "EBF1", "SPI1",
                "SMAD3", "DNMT1", "KLF2", "MYC"]


def load_panel():
    M = sp.load_npz(os.path.join(DATA_DIR, "panel_matrix.npz"))
    genes = [g for g in open(os.path.join(DATA_DIR, "panel_genes.txt"), encoding="utf-8").read().split("\n") if g]
    meta = pd.read_csv(os.path.join(DATA_DIR, "panel_meta.tsv"), sep="\t")
    return M, genes, meta


def build_regulons(genes, min_targets=5, max_targets=None):
    gi = {g: i for i, g in enumerate(genes)}
    # CollecTRI（OmniPath；第三方资源、不随仓库分发）：缺失时经 code/lib/collectri.py
    # 自动获取；文件已存在时读取路径与结果与原先完全一致。
    ctri = load_collectri(data_dir=DATA_DIR)
    ctri = ctri[["source_genesymbol", "target_genesymbol"]].dropna()
    reg = {}
    for tf, sub in ctri.groupby("source_genesymbol"):
        tg = sorted({t for t in sub["target_genesymbol"] if t in gi})
        if len(tg) >= min_targets and (max_targets is None or len(tg) <= max_targets):
            reg[tf] = tg
    return reg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-fig", action="store_true")
    ap.add_argument("--k-null", type=int, default=1000,
                    help="经验零分布：为每个调控子靶数生成 K 个同规模随机基因集")
    ap.add_argument("--no-partial", action="store_true",
                    help="跳过「剔除单元级全局协变量」的偏相关稳健性分析")
    a = ap.parse_args()

    print("=" * 74)
    print("step2 上游调控子推断与轴基因关联")
    print("=" * 74)
    M, genes, meta = load_panel()
    gi = {g: i for i, g in enumerate(genes)}

    # 细胞类型标签（step1 产物）
    asg = pd.read_csv(os.path.join(RES_DIR, "cell_type_assignment.tsv"), sep="\t")
    assert len(asg) == M.shape[1], "cell_type_assignment 与面板细胞数不一致，请先跑 step1"
    meta = meta.merge(asg[["cell_idx", "cell_type", "b_subset"]], on="cell_idx", how="left")
    print("细胞类型计数：")
    print(meta["cell_type"].value_counts().to_string())

    # ---- 限定 B 区室 ----
    bmask = meta["cell_type"].isin(["B", "Plasma"]).to_numpy()
    bmeta = meta.loc[bmask].reset_index(drop=True)
    Mb = M[:, np.where(bmask)[0]].tocsr()
    print(f"\nB 区室细胞数 = {Mb.shape[1]:,}")

    # ---- 基因 z 标准化（B 区室内） ----
    mu = np.asarray(Mb.mean(axis=1)).ravel()
    sq = np.asarray(Mb.multiply(Mb).mean(axis=1)).ravel()
    sd = np.sqrt(np.maximum(sq - mu ** 2, 1e-8))
    coo = Mb.tocoo()
    Z = sp.csr_matrix(((coo.data - mu[coo.row]) / sd[coo.row], (coo.row, coo.col)),
                      shape=Mb.shape, dtype=np.float32)
    np.clip(Z.data, -10, 10, out=Z.data)

    # ---- 调控子活性（ULM） ----
    reg = build_regulons(genes)
    print(f"可用调控子（面板内靶基因 ≥5、上限 不限）：{len(reg)}")
    tfs = sorted(reg)
    W = sp.lil_matrix((len(tfs), len(genes)), dtype=np.float32)
    for i, tf in enumerate(tfs):
        idx = [gi[t] for t in reg[tf]]
        W[i, idx] = 1.0 / len(idx)
    W = W.tocsr()
    ACT = W.dot(Z)                      # tfs x cells
    ACT = np.asarray(ACT.todense(), dtype=np.float32)
    act_df = pd.DataFrame(ACT.T, columns=tfs)
    print(f"调控子活性矩阵：{act_df.shape}")

    # ---- 轴基因/复合轴评分（B 区室，log1p CPM） ----
    ax_idx = [gi[g] for g in AXIS if g in gi]
    AX = np.asarray(Mb[ax_idx, :].todense(), dtype=np.float32).T   # cells x genes
    ax_names = [g for g in AXIS if g in gi]
    axda = pd.DataFrame(AX, columns=ax_names)
    axda["AxisComposite"] = axda[ax_names].mean(axis=1)

    # ---- 伪批量：按 (供者 × B 亚群) ----
    units = bmeta[["sample", "group", "b_subset"]].copy()
    assert (units["b_subset"] != "None").all(), "B 区室内存在未归入 B 亚群的细胞"
    act_df["_u"] = units["sample"].astype(str) + "|" + units["b_subset"].astype(str)
    axda["_u"] = act_df["_u"].values
    act_pb = act_df.groupby("_u").mean(numeric_only=True)
    ax_pb = axda.groupby("_u").mean(numeric_only=True)
    nunits = len(act_pb)
    print(f"伪批量单元数 = {nunits}")

    # ---- 关联检验：在每个 B 亚群内，跨供者 Spearman（复合轴评分 + 各单基因） ----
    rows = []
    for bs in SUBSETS:
        sel = [u for u in act_pb.index if u.endswith("|" + bs)]
        if len(sel) < 8:
            print(f"  [skip] {bs} 单元不足 ({len(sel)})")
            continue
        A = act_pb.loc[sel, tfs].to_numpy()
        Y = ax_pb.loc[sel].to_numpy()
        ynames = list(ax_pb.columns)
        for j, yn in enumerate(ynames):
            yv = Y[:, j]
            if np.nanstd(yv) == 0:
                continue
            for i, tf in enumerate(tfs):
                xv = A[:, i]
                if np.nanstd(xv) == 0:
                    continue
                rho, p = spearmanr(xv, yv, nan_policy="omit")
                if np.isnan(rho):
                    continue
                rows.append({"b_subset": bs, "n_donors": len(sel), "tf": tf,
                             "target": yn, "rho": float(rho), "p": float(p),
                             "n_targets": len(reg[tf])})
    res = pd.DataFrame(rows)
    res["q"] = np.nan
    for bs, sub in res.groupby("b_subset"):
        res.loc[sub.index, "q"] = multipletests(sub["p"].to_numpy(), method="fdr_bh")[1]

    # ---- 按「每个调控子自身靶数」匹配的经验零分布 ----
    # 做法：先把面板的 z 矩阵按伪批量单元预聚合为 D（基因 × 单元，等价于"先算每细胞
    #       活性再按单元取均值"），再对每个规模 s 抽 K 个 s 基因随机集，其活性 =
    #       D[随机集].mean(axis=0)，与复合轴评分跨单元 Spearman → 该规模的 ρ 零分布。
    #       这样每个调控子都只与"同样靶数"的随机集比较，量纲对齐。
    uid = act_pb.index.tolist()
    u_of = {u: i for i, u in enumerate(uid)}
    ux = np.array([u_of[u] for u in axda["_u"]], dtype=np.int32)
    C = sp.csr_matrix((np.ones(len(ux), dtype=np.float32), (ux, np.arange(len(ux)))),
                      shape=(len(uid), len(ux)))
    D = np.asarray((C.dot(Z.T)).toarray().T, dtype=np.float64)
    D = D / np.asarray(C.sum(axis=1)).ravel()[None, :]           # 基因 × 单元

    # 与"逐细胞活性"路径的一致性自检（数学等价；仅近并列项可能有浮点级差异）
    dev = 0.0
    for bs in SUBSETS:
        sel = [u for u in act_pb.index if u.endswith("|" + bs)]
        yv = ax_pb.loc[sel, "AxisComposite"].to_numpy()
        aix = np.array([u_of[u] for u in sel])
        for tf in tfs:
            r1 = spearmanr(act_pb.loc[sel, tf].to_numpy(), yv, nan_policy="omit")[0]
            r2 = spearmanr(D[[gi[t] for t in reg[tf]]][:, aix].mean(axis=0), yv,
                           nan_policy="omit")[0]
            if not (np.isnan(r1) or np.isnan(r2)):
                dev = max(dev, abs(r1 - r2))
    print(f"  [自检] 预聚合路径 vs 逐细胞路径  max|Δρ| = {dev:.3e}")

    rng = np.random.default_rng(0)
    K = a.k_null
    sizes = sorted({len(v) for v in reg.values()})
    print(f"  靶数档位数 = {len(sizes)}（{min(sizes)}–{max(sizes)}）；每档 K = {K}")

    def spearman_rows(A, yv):
        ra = rankdata(A, axis=1)
        ra = ra - ra.mean(axis=1, keepdims=True)
        ry = rankdata(yv)
        ry = ry - ry.mean()
        den = np.sqrt((ra ** 2).sum(axis=1)) * np.sqrt((ry ** 2).sum())
        with np.errstate(invalid="ignore", divide="ignore"):
            r = (ra @ ry) / den
        return r[np.isfinite(r)]          # 常数列（该随机集在单元间无变异）不计入零分布

    nulls = {}
    for bs in SUBSETS:
        sel = [u for u in act_pb.index if u.endswith("|" + bs)]
        if len(sel) < 8:
            continue
        yv = ax_pb.loc[sel, "AxisComposite"].to_numpy()
        aix = np.array([u_of[u] for u in sel])
        Dv = D[:, aix]
        per = {}
        for s in sizes:
            idx = np.array([rng.choice(len(genes), size=s, replace=False) for _ in range(K)])
            per[s] = spearman_rows(Dv[idx].mean(axis=1), yv)
        nulls[bs] = per
        med = np.median([np.median(per[s]) for s in sizes])
        nmin = min(len(per[s]) for s in sizes)
        print(f"  [按靶数匹配零分布] {bs}: 各档中位 ρ 的中位={med:.3f}；"
              f"有效随机集最少 {nmin}/{K}（常数列已剔除）")

    # ---- 汇总：每个 TF 与复合轴评分 ----
    comp = res[res["target"] == "AxisComposite"].copy()
    comp["null_med"] = np.nan
    comp["null_q95"] = np.nan
    comp["null_max"] = np.nan
    comp["emp_p"] = np.nan
    for bs, g in comp.groupby("b_subset"):
        for i, r in g.iterrows():
            nd = nulls[bs][int(r["n_targets"])]
            if len(nd) == 0:
                continue
            comp.loc[i, "null_med"] = float(np.median(nd))
            comp.loc[i, "null_q95"] = float(np.percentile(nd, 95))
            comp.loc[i, "null_max"] = float(nd.max())
            comp.loc[i, "emp_p"] = float((nd >= r["rho"]).mean())
    comp["q_emp"] = np.nan
    for bs, g in comp.groupby("b_subset"):
        comp.loc[g.index, "q_emp"] = multipletests(g["emp_p"].to_numpy(), method="fdr_bh")[1]
    comp = comp.sort_values(["b_subset", "rho"], ascending=[True, False])
    comp.to_csv(os.path.join(RES_DIR, "regulon_axis_composite.tsv"), sep="\t", index=False)
    res.to_csv(os.path.join(RES_DIR, "regulon_axis_correlation.tsv"), sep="\t", index=False)

    # ---- 稳健性控制：剔除「单元级全局协变量」后的偏相关 ----
    # 动机：面板基因在单元间共享一个很强的全局成分（文库复杂度 + 共同转录产出），
    #       它使任何基因集（含随机集）与轴评分都得到高 ρ。本节把该全局协变量
    #       g(u) = 面板全基因在单元 u 的平均 z 从**两边**同时线性剔除后再算 Spearman；
    #       零分布同样先对 g 取残差，保证与观测同口径。
    part_rows = []
    part = pd.DataFrame()
    if not a.no_partial:
        for bs in SUBSETS:
            sel = [u for u in act_pb.index if u.endswith("|" + bs)]
            if len(sel) < 8:
                continue
            aix = np.array([u_of[u] for u in sel])
            Dv = D[:, aix]                                  # 基因 × 单元
            G = Dv.mean(axis=0)                             # 单元级全局协变量
            Gc = G - G.mean()
            denG = float((Gc ** 2).sum())
            yv = ax_pb.loc[sel, "AxisComposite"].to_numpy()

            def _resid(v):
                v = np.asarray(v, dtype=np.float64)
                c = v - v.mean()
                if denG <= 0:
                    return c
                return c - (Gc @ c) / denG * Gc

            yr = _resid(yv)
            per = {}
            for s in sizes:
                idx = np.array([rng.choice(len(genes), size=s, replace=False) for _ in range(K)])
                act = Dv[idx].mean(axis=1)                  # K × 单元
                actc = act - act.mean(axis=1, keepdims=True)
                beta = (actc @ Gc) / denG if denG > 0 else np.zeros(len(actc))
                per[s] = spearman_rows(actc - beta[:, None] * Gc[None, :], yr)
            npairs = 0
            for tf in tfs:
                x = Dv[[gi[t] for t in reg[tf]]].mean(axis=0)
                xr = _resid(x)
                if np.std(xr) == 0:
                    continue
                rho_p = spearman_rows(xr[None, :], yr)
                if len(rho_p) == 0:
                    continue
                nd = per[len(reg[tf])]
                if len(nd) == 0:
                    continue
                part_rows.append({"b_subset": bs, "tf": tf, "n_targets": len(reg[tf]),
                                  "rho_partial": float(rho_p[0]),
                                  "emp_p": float((nd >= rho_p[0]).mean()),
                                  "null_med": float(np.median(nd)),
                                  "null_max": float(nd.max())})
                npairs += 1
            print(f"  [全局协变量控制] {bs}: 可算组合 {npairs}"
                  f"（随机集残差后中位 ρ ≈ {np.median(np.concatenate([per[s] for s in sizes])):.3f}）")
        part = pd.DataFrame(part_rows)
        if len(part):
            part["q_emp"] = np.nan
            for bs, g in part.groupby("b_subset"):
                part.loc[g.index, "q_emp"] = multipletests(g["emp_p"].to_numpy(), method="fdr_bh")[1]
            part = part.sort_values(["q_emp", "rho_partial"], ascending=[True, False])
            part.to_csv(os.path.join(RES_DIR, "regulon_axis_partial.tsv"), sep="\t", index=False)
            print(f"  偏相关：q_emp < 0.10 共 {int((part['q_emp'] < 0.10).sum())} / {len(part)}"
                  f"（靶数上限 {int(part.loc[part['q_emp'] < 0.10, 'n_targets'].max()) if (part['q_emp'] < 0.10).any() else 0}）")

    n_combos = len(comp)
    n_emp05 = int((comp["emp_p"] < 0.05).sum())
    n_q05 = int((comp["q_emp"] < 0.05).sum())
    n_q10 = int((comp["q_emp"] < 0.10).sum())
    print("\n复合轴评分关联 top（按 B 亚群；emp_p = 同靶数随机基因集中 ≥ 观测 ρ 的比例）：")
    for bs, sub in comp.groupby("b_subset"):
        print(f"  [{bs}]")
        print(sub.head(6)[["tf", "rho", "emp_p", "q_emp", "n_targets"]].to_string(index=False))
    print(f"\n组合数 = {n_combos}；emp_p < 0.05 共 {n_emp05} 个（随机期望 {0.05*n_combos:.1f} 个）")
    print(f"BH-FDR：q_emp < 0.05 共 {n_q05} 个；q_emp < 0.10 共 {n_q10} 个")

    # 预先指定的候选 TF 集（描述性，不改变上述判定）
    cand = comp[comp["tf"].isin(CYTOKINE_TFS)]
    print(f"预先指定候选集：{len(cand)} 个组合中 emp_p<0.05 的 = {int((cand['emp_p']<0.05).sum())}"
          f"（期望 {0.05*len(cand):.2f}）")
    if (cand["emp_p"] < 0.05).any():
        print(cand[cand["emp_p"] < 0.05][["tf", "b_subset", "n_targets", "rho", "emp_p"]]
              .to_string(index=False))

    # 细胞因子 TF 专项
    cyto = res[res["tf"].isin(CYTOKINE_TFS)]
    cyto.to_csv(os.path.join(RES_DIR, "regulon_cytokine_tfs.tsv"), sep="\t", index=False)

    # TF 活性在 B 亚群间差异（描述性）
    sub_mean = pd.DataFrame({
        bs: act_pb.loc[[u for u in act_pb.index if u.endswith("|" + bs)]].mean()
        for bs in SUBSETS})
    sub_mean.to_csv(os.path.join(RES_DIR, "regulon_activity_by_bsubset.tsv"), sep="\t")

    if not a.no_fig:
        make_fig(res, comp, act_pb, ax_pb, tfs, nulls, sizes)

    summary = {"n_b_cells": int(Mb.shape[1]), "n_regulons": len(tfs),
               "target_rule": {"min_targets": 5, "max_targets": None,
                               "note": "2026-10-01 起不设上限；按自身靶数匹配的零分布"},
               "target_median": int(np.median([len(v) for v in reg.values()])),
               "target_max": int(max(len(v) for v in reg.values())),
               "n_combos": int(n_combos), "n_passed_emp05": n_emp05,
               "n_passed_q05": n_q05, "n_passed_q10": n_q10,
               "k_null": int(K), "n_units": int(nunits),
               "b_subset_counts": bmeta["b_subset"].value_counts().to_dict(),
               "preaggregation_max_dev": float(dev),
               "partial_control": {
                   "note": "剔除单元级全局协变量（面板全基因 z 的单元均值）后的偏相关",
                   "n_pairs": int(len(part)),
                   "n_passed_q10": int((part["q_emp"] < 0.10).sum()) if len(part) else 0,
                   "max_targets_passing_q10": (int(part.loc[part["q_emp"] < 0.10, "n_targets"].max())
                                               if len(part) and (part["q_emp"] < 0.10).any() else 0),
                   "file": "regulon_axis_partial.tsv"},
               "top_composite": comp.groupby("b_subset").head(5)[
                   ["b_subset", "tf", "rho", "emp_p", "q_emp", "n_targets"]].to_dict("records")}
    with open(os.path.join(RES_DIR, "step2_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print("\n[done] step2")


def make_fig(res, comp, act_pb, ax_pb, tfs, nulls, sizes):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(15, 5.0))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.0, 1.0], wspace=0.42)

    # (a) 观测 ρ 与「同靶数随机集」零分布包络：看调控子是否真的超出随机期望
    ax = fig.add_subplot(gs[0, 0])
    bs_show = "B_naive" if "B_naive" in nulls else list(nulls)[0]
    per = nulls[bs_show]
    xs = sorted(sizes)
    med = [np.median(per[s]) for s in xs]
    q95 = [np.percentile(per[s], 95) for s in xs]
    ax.plot(xs, med, color="#7f7f7f", lw=1.2, label="random sets, median ρ")
    ax.plot(xs, q95, color="#7f7f7f", lw=1.2, ls="--", label="random sets, 95th pct")
    ax.fill_between(xs, med, q95, color="#d9d9d9", alpha=0.55, zorder=0)
    d = comp[comp["b_subset"] == bs_show]
    ax.scatter(d["n_targets"], d["rho"], s=7, color="#c0392b", alpha=0.75,
               label="observed regulons", zorder=3)
    ax.set_xscale("log")
    ax.set_xlabel("regulon size (targets in panel)")
    ax.set_ylabel("Spearman ρ vs composite axis score")
    ax.set_title("Regulon activity vs axis score (%s)\n"
                 "target-number-matched random-set envelope" % bs_show.replace("_", " "),
                 fontsize=10.5, fontweight="bold")
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    ax.grid(alpha=0.25)

    # (b) 预先指定候选 TF 热图（naive）
    ax = fig.add_subplot(gs[0, 1])
    show = [t for t in CYTOKINE_TFS if t in tfs]
    m = comp[comp["b_subset"] == "B_naive"].set_index("tf").reindex(show)["rho"]
    m = m.dropna()
    if len(m):
        im = ax.imshow(m.values.reshape(-1, 1), cmap="YlOrRd", vmin=0.4, vmax=0.95, aspect="auto")
        ax.set_yticks(range(len(m))); ax.set_yticklabels(m.index, fontsize=8)
        ax.set_xticks([0]); ax.set_xticklabels(["axis\n(naive B)"], fontsize=8)
        for k, v in enumerate(m.values):
            ax.text(0, k, f"{v:.2f}", ha="center", va="center", fontsize=6.5,
                    color="white" if v > 0.78 else "black")
        fig.colorbar(im, ax=ax, fraction=0.06, pad=0.03).set_label("Spearman ρ", fontsize=7.5)
    ax.set_title("Pre-specified candidate TFs\n(cytokine / B-cell / epigenetic machinery)",
                 fontsize=10.5, fontweight="bold")

    # (c) 候选 TF 中 ρ 最高者 vs 轴评分散点（naive）
    ax = fig.add_subplot(gs[0, 2])
    sel = [u for u in act_pb.index if u.endswith("|B_naive")]
    cand = list(comp[comp["b_subset"] == "B_naive"].head(3)["tf"]) if not comp[comp["b_subset"] == "B_naive"].empty else []
    cand = [c for c in cand if c in act_pb.columns]
    yv = ax_pb.loc[sel, "AxisComposite"].to_numpy()
    for i, tf in enumerate(cand):
        xv = act_pb.loc[sel, tf].to_numpy()
        ax.scatter(xv, yv, s=18, alpha=0.8, label=f"{tf} (ρ={spearmanr(xv, yv)[0]:.2f})")
    ax.set_xlabel("regulon activity"); ax.set_ylabel("composite axis score")
    ax.legend(fontsize=7, frameon=False)
    ax.set_title("Highest-ρ regulators vs axis\n(naive B, per-donor)", fontsize=10.5, fontweight="bold")
    ax.grid(alpha=0.3)

    fig.tight_layout()
    for ext in ("png", "pdf"):
        out = os.path.join(FIG_DIR, f"Figure_S6_scRNA_regulons.{ext}")
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print(f"  图已写：{out}")
    try:
        from PIL import Image
        im = Image.open(os.path.join(FIG_DIR, "Figure_S6_scRNA_regulons.png"))
        im.save(os.path.join(FIG_DIR, "Figure_S6_scRNA_regulons.tiff"),
                format="TIFF", compression="tiff_lzw", dpi=(600, 600))
        print("  TIFF 已写")
    except Exception as e:
        print(f"  [警告] TIFF 失败：{e}")


if __name__ == "__main__":
    main()
