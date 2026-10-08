# -*- coding: utf-8 -*-
"""code/_paths.py —— 复现包统一路径解析（"开箱即用"）

目的
----
`code/` 下的脚本在**任意机器、任意克隆目录**下无需改动即可运行：所有脚本
通过本模块取路径，不再硬编码任何本机绝对路径。

目录约定
--------
REPO        复现包根目录（含 `code/`、`data/`、`results/`）
  ├ CODE    `<REPO>/code`
  ├ DATA    `<REPO>/data`          —— 各证据层的处理后结果表（随仓库分发）
  ├ RESULTS `<REPO>/results`       —— 分析产出（图表/汇总，随仓库分发）
  └ LIB     `<REPO>/code/lib`      —— 随包分发的辅助模块（如 tabix_remote.py）
RAWDATA     外部原始数据根（原分析的 `00_rawdata/`）—— **不随仓库分发**，
            需按论文 Data Availability Statement 所列 accession 自行下载到此。
SUBMISSION  投稿材料目录（原分析的 `阶段4.5_复现包/submission/`）—— 属投稿件，
            **不随仓库分发**；仅 `code/review_cleanup_*` 的审计脚本会用到。

环境变量（均可选，优先级最高）
------------------------------
  IGAN_REPRO_ROOT       复现包根目录
  IGAN_RAWDATA_ROOT     外部原始数据根（原 `00_rawdata/`）
  IGAN_SUBMISSION_DIR   投稿材料目录（原 `阶段4.5_复现包/submission/`）

RAWDATA 的默认解析顺序（取第一个真实存在的目录）：
  $IGAN_RAWDATA_ROOT → <REPO>/00_rawdata → <REPO>/rawdata
                     → <REPO>/../00_rawdata → <REPO>/../rawdata
都找不到时仍返回 <REPO>/../00_rawdata，使调用方能给出明确的"请先下载原始数据"提示。

用法
----
    import _paths
    df = pd.read_csv(_paths.data("M2_lead_eqtl", "lead_eQTL_matrix.tsv"))
    raw = _paths.raw("eQTL_Catalogue", "region_cache")          # 外部原始数据
    out = _paths.results("figures", "Fig1_v4.2.png")
    _paths.check()          # 自检：打印各根目录及存在性

历史脚本迁移
------------
早期脚本使用原工程树的绝对路径（如 `…\\IgAshenbing\\阶段2\\主证据1_lead_eQTL`）。
`_paths.LEGACY(...)` 与 `legacy(...)` 保留这些旧写法并映射到上表位置，
既可调用（`LEGACY("00_rawdata", "GEO_M1")`）又是一个字符串
（`LEGACY` == `<REPO>/..`，可直接赋给 `subprocess(cwd=...)`）。
"""

import os
import re

__all__ = [
    "REPO", "CODE", "DATA", "RESULTS", "LIB", "RAWDATA", "SUBMISSION",
    "repo", "code", "data", "results", "lib", "raw", "submission",
    "legacy", "LEGACY", "at", "check",
]


def _abspath(p):
    return os.path.abspath(os.path.expanduser(p)) if p else None


def _first_dir(*cands):
    for c in cands:
        if c and os.path.isdir(c):
            return c
    return None


# --------------------------------------------------------------------------
# 1. 复现包根目录：本文件即 <REPO>/code/_paths.py
# --------------------------------------------------------------------------
REPO = _abspath(os.environ.get("IGAN_REPRO_ROOT")) or \
       os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

CODE    = os.path.join(REPO, "code")
DATA    = os.path.join(REPO, "data")
RESULTS = os.path.join(REPO, "results")
LIB     = os.path.join(CODE, "lib")

# 原工程树根（`阶段1/…`、`00_rawdata/` 所在目录）；仅用于旧路径映射与审计脚本
_ENGINE = os.path.dirname(REPO)             # 本包在原始工程中的父目录（通常 = …/阶段4.5_复现包）
_PROJECT = os.path.dirname(_ENGINE)         # 原始工程根（通常 = …/IgAshenbing）

# --------------------------------------------------------------------------
# 2. 外部原始数据根（不随仓库分发）
# --------------------------------------------------------------------------
RAWDATA = _first_dir(
    _abspath(os.environ.get("IGAN_RAWDATA_ROOT")),
    os.path.join(REPO, "00_rawdata"),
    os.path.join(REPO, "rawdata"),
    os.path.join(_ENGINE, "00_rawdata"),
    os.path.join(_ENGINE, "rawdata"),
    os.path.join(_PROJECT, "00_rawdata"),
    os.path.join(_PROJECT, "rawdata"),
) or os.path.join(_ENGINE, "00_rawdata")

# --------------------------------------------------------------------------
# 3. 投稿材料目录（不随仓库分发）
# --------------------------------------------------------------------------
SUBMISSION = _first_dir(
    _abspath(os.environ.get("IGAN_SUBMISSION_DIR")),
    os.path.join(_ENGINE, "submission"),
    os.path.join(_ENGINE, "阶段4.5_复现包", "submission"),
    os.path.join(REPO, "submission"),
) or os.path.join(_ENGINE, "阶段4.5_复现包", "submission")


# --------------------------------------------------------------------------
# 4. 取路径
# --------------------------------------------------------------------------
def repo(*p):
    return os.path.join(REPO, *p) if p else REPO


def code(*p):
    return os.path.join(CODE, *p) if p else CODE


def data(*p):
    return os.path.join(DATA, *p) if p else DATA


def results(*p):
    return os.path.join(RESULTS, *p) if p else RESULTS


def lib(*p):
    return os.path.join(LIB, *p) if p else LIB


def raw(*p):
    return os.path.join(RAWDATA, *p) if p else RAWDATA


def submission(*p):
    return os.path.join(SUBMISSION, *p) if p else SUBMISSION


# --------------------------------------------------------------------------
# 5. 旧工程路径 → 复现包内位置（历史脚本迁移用）
# --------------------------------------------------------------------------
# 长键优先匹配（`阶段4.5_复现包/reproducibility` 先于 `阶段4.5_复现包`）
_TABLE = {
    "00_rawdata":                        raw,
    "阶段1/M1_GEO表达锚点":               lambda *s: data("M1_geo", *s),
    "阶段1/M2a_eQTL语境":                 lambda *s: data("M2a_eqtl", *s),
    "阶段2/主证据1_lead_eQTL":            lambda *s: data("M2_lead_eqtl", *s),
    "阶段2/主证据2_共享结构":             lambda *s: data("M2b_susie_coloc", *s),
    "阶段2/主证据3_机制轴不驱动疾病易感": lambda *s: data("M2c_disease", *s),
    "阶段3/M4_中介":                      lambda *s: data("M4_mediation", *s),
    "阶段3/M5_东亚":                      lambda *s: data("M5_eastasia", *s),
    "阶段3/M6_干预":                      lambda *s: data("M6_intervention", *s),
    "阶段4/主图/scripts":                 lambda *s: code("fig_scripts", *s),
    "阶段4/主图/out":                     lambda *s: results("figures", *s),
    "阶段4.5_复现包/reproducibility":      repo,
    "阶段4.5_复现包/submission":           submission,
    "阶段4.5_复现包":                      lambda *s: os.path.join(_ENGINE, *s) if s else _ENGINE,
    "submission":                         submission,
    "reproducibility":                    repo,
    ".workbuddy/tmp/tabix":               lib,
}
# 未随仓库分发的原工程阶段目录（`阶段4/初稿`、`阶段4/主图/data_check` 等）：
# 只在原始工程树内可解析，供审计脚本使用。
for _st in ("阶段1", "阶段2", "阶段3", "阶段4"):
    _TABLE[_st] = (lambda st: (lambda *s: os.path.join(_PROJECT, st, *s)))(_st)
_KEYS = sorted(_TABLE, key=len, reverse=True)


def legacy(*parts):
    """把原工程树相对路径（`\\` 或 `/` 分隔）映射到复现包内等效位置。

    >>> legacy(r"阶段2\\主证据1_lead_eQTL", "lead_eQTL_matrix.tsv")
    '<REPO>/data/M2_lead_eqtl/lead_eQTL_matrix.tsv'
    >>> legacy("00_rawdata", "GEO_M1")          # 外部原始数据 → RAWDATA
    """
    raw_join = "/".join(str(p) for p in parts)
    segs = [s for s in re.split(r"[\\/]+", raw_join) if s]
    for k in _KEYS:
        ks = k.split("/")
        if segs[:len(ks)] == ks:
            return _TABLE[k](*segs[len(ks):])
    # 未收录的旧路径（如投稿件 `阶段4/初稿`）：回落到原工程树根下
    return os.path.join(_ENGINE, *segs) if segs else _ENGINE


class _LegacyCallable(str):
    """既是字符串（可直接当路径 / `cwd=` 用），又可调用（等价 `legacy(前缀…, …)`）。"""

    def __new__(cls, value, prefix=()):
        obj = super().__new__(cls, value)
        obj._prefix = tuple(prefix)
        return obj

    def __call__(self, *parts):
        return legacy(*self._prefix, *parts)


LEGACY = _LegacyCallable(_ENGINE)          # 字符串值 = 原工程树根


def at(*prefix):
    """返回"已绑定前缀"的 LEGACY：`at("submission", "manuscript")("x.md")`
    等价 `legacy("submission", "manuscript", "x.md")`，且其字符串值即映射后的目录。"""
    return _LegacyCallable(legacy(*prefix), prefix)


# --------------------------------------------------------------------------
# 6. 自检
# --------------------------------------------------------------------------
def check(verbose=True):
    """打印各根目录与关键输入的存在性，返回 0 问题数。"""
    problems = 0
    rows = [
        ("REPO", REPO, True),
        ("CODE", CODE, True),
        ("DATA", DATA, True),
        ("RESULTS", RESULTS, True),
        ("LIB", LIB, True),
        ("RAWDATA (外部原始数据，需自行下载)", RAWDATA, False),
        ("SUBMISSION (投稿件，不随仓库分发)", SUBMISSION, False),
    ]
    for label, p, must in rows:
        ok = os.path.isdir(p)
        flag = "OK " if ok else ("MISS" if must else "--  ")
        if must and not ok:
            problems += 1
        if verbose:
            print(f"[{flag}] {label:<38} {p}")
    for name in ("M2a_eqtl", "M1_geo", "M2_lead_eqtl", "M2b_susie_coloc",
                 "M2c_disease", "M4_mediation", "M5_eastasia",
                 "M6_intervention", "M7_ancestry_increment", "sc_regulatory"):
        p = os.path.join(DATA, name)
        if not os.path.isdir(p):
            problems += 1
            if verbose:
                print(f"[MISS] data/{name}")
    if verbose:
        print(f"-- {problems} problem(s)")
    return problems


if __name__ == "__main__":
    raise SystemExit(check())
