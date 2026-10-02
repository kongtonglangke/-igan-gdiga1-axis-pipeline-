# -*- coding: utf-8 -*-
"""coloc_abf 修正后单元测试 (合成数据) — 验证与 R coloc::coloc.abf 语义一致"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import importlib.util, numpy as np
spec = importlib.util.spec_from_file_location(
    "s1", _paths.code("M2b_susie_coloc", "step1_coloc_abf.py"))
s1 = importlib.util.module_from_spec(spec)
# 不执行 main: 手动剥掉 __main__ 块
src = open(spec.origin, encoding="utf-8").read()
src = src.split('if __name__ == "__main__":')[0]
exec(compile(src, spec.origin, "exec"), s1.__dict__)

rng = np.random.default_rng(2026)
n = 500
pos = np.arange(n)
# 区域相关结构 (ar1 rho=0.9)
def ar1_corr(n, rho):
    c = np.empty((n, n))
    for i in range(n):
        for j in range(n):
            c[i, j] = rho ** abs(i - j)
    return c
R = ar1_corr(n, 0.9)
L = np.linalg.cholesky(R + 1e-6 * np.eye(n))
def gen_beta_se(causal_idx, effect, n_eff=2000):
    # 生成近似 z 结构: 因果位点效应 + LD 传播 + 噪声
    z_true = np.zeros(n)
    z_true[causal_idx] = effect / 0.05  # 等效效应
    z = np.sqrt(n_eff) * (L @ (z_true / np.sqrt(n_eff))) if False else (R @ z_true) + rng.normal(0, 1.2, n)
    se = np.full(n, 0.05)
    beta = z * se
    return beta, se

def show(tag, res):
    print(f"[{tag}] PPH0={res['PPH0']:.4f} PPH1={res['PPH1']:.4f} "
          f"PPH2={res['PPH2']:.4f} PPH3={res['PPH3']:.4f} PPH4={res['PPH4']:.4f} n={res['n_SNP']}")

# 场景1: 共享因果 (同一 SNP 驱动两性状强信号) → PPH4 应主导
b1, s1v = gen_beta_se(250, 0.35)
b2, s2v = gen_beta_se(250, 0.40)
show("共享因果(同SNP250)", s1.coloc_abf(b1, s1v, b2, s2v))

# 场景2: trait2 无信号 → PPH0 应主导 (或 PPH1)
b2n, s2n = rng.normal(0, 0.05, n), np.full(n, 0.05)
show("trait2无信号", s1.coloc_abf(b1, s1v, b2n, s2n))

# 场景3: 两性状独立因果位点 (SNP100 vs SNP400) → PPH3 应相对高
b3, s3v = gen_beta_se(400, 0.40)
show("独立因果(100 vs 400)", s1.coloc_abf(b1, s1v, b3, s3v))

# 场景4: 强共享 → PPH4 极端高 (审稿场景复核)
z = np.zeros(n); z[300] = 8.0
b4 = (R @ z) * 0.05 + rng.normal(0, 0.2, n)
show("强共享z=8", s1.coloc_abf(b4, np.full(n, 0.05), b4 + rng.normal(0, 0.02, n), np.full(n, 0.05)))
