#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
step4_supp_tables.py — 汇总 step2/step3 产出为投稿用补充表（Table S9 / Table S10）
================================================================================
把分析树里的冻结结果导出为投稿材料中的两张补充表，避免手工拼表：

  Table S9  ← results/sc_regulatory/TableS9_methylation_sensitive_TFs.tsv      （逐字节复制）
  Table S10 ← results/sc_regulatory/regulon_axis_correlation.tsv（逐基因行）
              + regulon_axis_composite.tsv（复合轴评分行的零分布与经验 P / q）

Table S10 会为 `target == AxisComposite` 的行追加 5 个列：
  axis_composite_null_median / _null_q95 / _null_max / _emp_p / _q_emp
逐基因行这几列留空（按自身靶数匹配的经验零分布只适用于复合轴评分检验）。

旧文件在覆盖前归档为 `<name>.before_<tag>_<date>.tsv`（留档不删）。

用法：
  python code/sc_regulatory/step4_supp_tables.py            # 导出
  python code/sc_regulatory/step4_supp_tables.py --check    # 只比对，不写盘
"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402
# --------------------------------------------------------------------------
import csv
import shutil
import sys

RES = _paths.results("sc_regulatory")
ADD = _paths.submission("additional_file_1")
TAG = "planB"; DATE = "20261001"

NULL_COLS = ["axis_composite_null_median", "axis_composite_null_q95",
             "axis_composite_null_max", "axis_composite_emp_p", "axis_composite_q_emp"]


def read_tsv(p):
    with open(p, encoding="utf-8", newline="") as f:
        rdr = csv.DictReader(f, delimiter="\t")
        return list(rdr), list(rdr.fieldnames)


def write_tsv(p, rows, cols):
    with open(p, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def build_s10():
    corr, ccols = read_tsv(_os.path.join(RES, "regulon_axis_correlation.tsv"))
    comp, _ = read_tsv(_os.path.join(RES, "regulon_axis_composite.tsv"))
    assert all(r["target"] == "AxisComposite" for r in comp), "composite 表含非 AxisComposite 行"
    smap = {}
    for r in comp:
        k = (r["b_subset"], r["tf"])
        assert k not in smap, f"composite 表重复键 {k}"
        smap[k] = (r["null_med"], r["null_q95"], r["null_max"], r["emp_p"], r["q_emp"])

    n_comp = sum(1 for r in corr if r["target"] == "AxisComposite")
    assert n_comp == len(comp), f"组合数不一致：corr {n_comp} vs comp {len(comp)}"
    assert {(r["b_subset"], r["tf"]) for r in corr if r["target"] == "AxisComposite"} == set(smap), \
        "键集合不一致"

    out_cols = ccols + NULL_COLS
    for r in corr:
        if r["target"] == "AxisComposite":
            vals = smap[(r["b_subset"], r["tf"])]
            for c, v in zip(NULL_COLS, vals):
                r[c] = v
        else:
            for c in NULL_COLS:
                r[c] = ""
    return corr, out_cols


def main():
    check = "--check" in sys.argv

    # ---------- Table S9：逐字节复制 ----------
    s9_src = _os.path.join(RES, "TableS9_methylation_sensitive_TFs.tsv")
    s9_dst = _os.path.join(ADD, "Table_S9_methylation_sensitive_TFs.tsv")
    same9 = (open(s9_src, "rb").read() == open(s9_dst, "rb").read()) \
        if _os.path.exists(s9_dst) else False
    print(f"Table S9 : {'一致' if same9 else '需更新'}  ({_os.path.getsize(s9_src):,} B)")

    # ---------- Table S10 ----------
    rows, cols = build_s10()
    s10_dst = _os.path.join(ADD, "Table_S10_regulon_axis_associations.tsv")
    same10 = False
    if _os.path.exists(s10_dst):
        old, ocols = read_tsv(s10_dst)
        same10 = (len(old) == len(rows) and ocols == cols
                  and all(o == n for o, n in zip(old, rows)))
    fill = sum(1 for r in rows if r["target"] == "AxisComposite")
    print(f"Table S10: {'一致' if same10 else '需更新'}  "
          f"({len(rows):,} 行 × {len(cols)} 列；AxisComposite 填充 {fill:,} 行)")

    if check:
        print("CHECK:", "ALL MATCH" if (same9 and same10) else "NEED UPDATE")
        return 0 if (same9 and same10) else 1

    if not same9:
        shutil.copyfile(s9_src, s9_dst)
        print(f"  已写 {_os.path.basename(s9_dst)}")
    if not same10:
        if _os.path.exists(s10_dst):
            arc = s10_dst.replace(".tsv", f".before_{TAG}_{DATE}.tsv")
            if not _os.path.exists(arc):
                shutil.copy2(s10_dst, arc)
                print(f"  旧版归档 → {_os.path.basename(arc)}")
        write_tsv(s10_dst, rows, cols)
        print(f"  已写 {_os.path.basename(s10_dst)}")

    # ---------- 头条数值自检（与正文一致） ----------
    key = {(r["b_subset"], r["tf"]): r for r in rows if r["target"] == "AxisComposite"}
    exp = {("B_naive", "ERF"): "0.0", ("B_naive", "STAT3"): "0.5", ("B_naive", "SP1"): "0.393"}
    for k, v in exp.items():
        got = key[k]["axis_composite_emp_p"]
        flag = "OK " if got == v else "!! "
        print(f"  [{flag}] {k[0]} / {k[1]}  empirical P = {got}（期望 {v}）")
    print("[done] step4")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
