#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
step1_immunexut_crossref.py — C1：东亚（日本）细胞类型 eQTL 交叉参照
================================================================================
目的：现稿 Limitations #5 指出"解读东亚信号的 eQTL 面板（OneK1K、GTEx）均为欧洲
      来源"。本步引入 **ImmuNexUT**（Ota et al. Cell 2021;184:3006-3021.e17）——
      416 名日本受试者、28 种免疫细胞亚型（含 Naive_B / USM_B / SM_B / DN_B /
      Plasmablast）的条件独立 eQTL（FDR<0.05）——作为**东亚细胞类型分辨**的
      交叉参照，检验五个轴基因是否在东亚 B 系亚型中为 eGene、以及 Gd-IgA1 lead
      变异是否出现在其条件独立信号中。

数据：00_rawdata/ImmuNexUT/E-GEAD-398.processed.zip（NBDC Human DB, E-GEAD-398）
      └─ 内含单个成员 conditional_eQTL_FDR0.05.tar（解压后约 4.1 GB）
      本脚本**流式**读取（zipfile → tarfile 'r|'），不落地 4.1 GB tar。

用法：
  # 先探针：列出 tar 成员名 + 每个成员首行（了解列结构），不解析
  python step1_immunexut_crossref.py --peek

  # 正式抽取：过滤五个轴基因的全部行 → 结果目录
  python step1_immunexut_crossref.py

输出（results/C1_eas_eqtl/）：
  immunexut_members.txt         tar 成员清单
  axis_hits_raw.tsv             轴基因命中行（member + 原始行）
  step1_summary.json            汇总（每亚型的轴基因 eGene、最小 P、lead 变异命中）
  TableS12_EAS_eQTL_crossref.tsv  （由 step2 汇总；本脚本只产出 raw + summary）
"""

# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import io
import os
import sys
import json
import tarfile
import zipfile

BASE = _paths.LEGACY
ZIP = BASE("00_rawdata", "ImmuNexUT", "E-GEAD-398.processed.zip")
RES = BASE("阶段4.5_复现包", "reproducibility", "results", "C1_eas_eqtl")
os.makedirs(RES, exist_ok=True)

# 五个轴基因（Ensembl 经 Ensembl REST 核对，2026-09-10）
AXIS = {
    "ENSG00000106392": "C1GALT1",
    "ENSG00000171155": "C1GALT1C1",
    "ENSG00000143641": "GALNT2",
    "ENSG00000119514": "GALNT12",
    "ENSG00000070731": "ST6GALNAC2",
}
SYMS = set(AXIS.values())
# 四个已发表 Gd-IgA1 lead 变异（含 rsID 与 hg19/hg38 位置两种匹配方式）
LEADS = ["rs13226913", "rs10238682", "rs7856182", "rs5910940"]

PEEK = "--peek" in sys.argv


def iter_tar_members():
    """流式产出 (member_name, fileobj) —— 不落地 tar。"""
    zf = zipfile.ZipFile(ZIP)
    tar_members = [n for n in zf.namelist() if n.lower().endswith(".tar")]
    if not tar_members:
        raise RuntimeError(f"zip 内未找到 .tar 成员；实际成员：{zf.namelist()}")
    member = tar_members[0]
    print(f"[zip] 成员 {member}（压缩 {zf.getinfo(member).compress_size/1e6:.1f} MB "
          f"→ 原始 {zf.getinfo(member).file_size/1e9:.2f} GB）", flush=True)
    raw = zf.open(member)
    tf = tarfile.open(fileobj=raw, mode="r|")   # 流式，不支持随机访问
    for ti in tf:
        if not ti.isfile():
            continue
        yield ti.name, tf.extractfile(ti)


def main():
    print("=" * 76)
    print("C1：ImmuNexUT（日本，28 亚型）轴基因 eQTL 交叉参照")
    print("=" * 76)

    names, sample_lines, n_members, n_bytes = [], {}, 0, 0
    hits_path = os.path.join(RES, "axis_hits_raw.tsv")
    fh_out = None if PEEK else open(hits_path, "w", encoding="utf-8")
    if fh_out:
        fh_out.write("member\tline\n")

    per_member_hits = {}
    lead_hits = []

    for name, fobj in iter_tar_members():
        n_members += 1
        names.append(name)
        if n_members % 25 == 0:
            print(f"  ... {n_members} members, {n_bytes/1e9:.2f} GB streamed", flush=True)
        kept = 0
        first = None
        for raw in fobj:
            n_bytes += len(raw)
            line = raw.decode("utf-8", "replace").rstrip("\n")
            if first is None:
                first = line
            up = line
            # 命中：任一轴基因 ENSG 或基因符号；或任一 lead rsID
            is_axis = any(e in up for e in AXIS) or any(s in up for s in SYMS)
            is_lead = False
            if not is_axis:
                for rs in LEADS:
                    if rs in up:
                        is_lead = True
                        lead_hits.append((name, line[:400]))
                        break
            if is_axis:
                kept += 1
                if fh_out:
                    fh_out.write(f"{name}\t{line}\n")
        if first is not None:
            sample_lines[name] = first
        if kept:
            per_member_hits[name] = kept
            print(f"  [hit] {name}: {kept} 行", flush=True)

    if fh_out:
        fh_out.close()

    with open(os.path.join(RES, "immunexut_members.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(names))
    if PEEK:
        peek = os.path.join(RES, "immunexut_peek.tsv")
        with open(peek, "w", encoding="utf-8") as f:
            for n in names[:60]:
                f.write(f"{n}\t{sample_lines.get(n,'')}\n")
        print(f"\n[peek] members={n_members} 表头样例 → {peek}")
        for n in names[:12]:
            print(f"  {n}  ||  {sample_lines.get(n,'')[:160]}")

    summary = {
        "n_members": n_members,
        "bytes_streamed": n_bytes,
        "n_axis_hit_lines": sum(per_member_hits.values()),
        "per_member_hit_counts": per_member_hits,
        "n_lead_hit_lines": len(lead_hits),
        "lead_hits": [{"member": m, "line": l} for m, l in lead_hits[:50]],
        "axis_ensg": AXIS,
    }
    with open(os.path.join(RES, "step1_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print("\n[done] C1 step1")
    print(f"  members={n_members}  streamed={n_bytes/1e9:.2f} GB  "
          f"axis hit lines={sum(per_member_hits.values())}  lead hit lines={len(lead_hits)}")


if __name__ == "__main__":
    main()
