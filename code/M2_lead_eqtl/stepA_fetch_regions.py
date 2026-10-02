# -*- coding: utf-8 -*-
"""阶段2 主证据1 —— 步骤A：远程区域抓取落盘缓存
4 lead Gd-IgA1 位点 × 5 个 eQTL 语境（OneK1K B_naive/B_memory/B_intermediate +
GTEx_v10 WholeBlood/KidneyCortex），按位点 ±500kb 区域远程拉取 all.tsv.gz 原始行，
落盘到 00_rawdata/eQTL_Catalogue/region_cache/{qtd}.{chrom}.tsv 供步骤B解析。
chr7 上 rs13226913(7207215) 与 rs10238682(7215386) 相距仅 8kb → 合并为同一区域一次抓取。
"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import os, sys, time, gzip, io, urllib.request, concurrent.futures as cf

sys.path.insert(0, _paths.LIB)   # 随包分发的 tabix_remote.py
from tabix_remote import parse_tbi, region_file_ranges, decompress_bgzf

TBI_DIR = _paths.raw("eQTL_Catalogue", "tbi")   # .tbi 索引缓存（缺失时自动下载）
CACHE_DIR = _paths.raw("eQTL_Catalogue", "region_cache")
os.makedirs(CACHE_DIR, exist_ok=True)

def tbi_path(qtd, qts):
    """{qtd}.tbi 索引路径；缺失时自动从 eQTL Catalogue 下载（第三方索引，不随仓库分发）。"""
    os.makedirs(TBI_DIR, exist_ok=True)
    p = os.path.join(TBI_DIR, f"{qtd}.tbi")
    if not os.path.exists(p):
        url = f"{BASE}/{qts}/{qtd}/{qtd}.all.tsv.gz.tbi"
        print(f"[tbi] 下载索引 {url}", flush=True)
        urllib.request.urlretrieve(url, p)
    return p


BASE = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats"
# (label, QTD, QTS)
DATASETS = [
    ("OneK1K_B_naive",        "QTD000608", "QTS000038"),
    ("OneK1K_B_memory",       "QTD000607", "QTS000038"),
    ("OneK1K_B_intermediate", "QTD000606", "QTS000038"),
    ("GTEx_v10_blood",        "QTD000356", "QTS000015"),
    ("GTEx_v10_kidney_cortex","QTD000261", "QTS000015"),
]
# (chrom, beg, end) 1-based 闭区间 —— 由 4 lead ±500kb 合并去重
REGIONS = [
    ("7", 6707215, 7715386),   # rs13226913 ±500kb ∪ rs10238682 ±500kb
    ("9", 98371030, 99371030), # rs7856182 ±500kb (GALNT12)
    ("X", 120124955, 121124955),  # rs5910940 ±500kb (C1GALT1C1)
]

def file_size(url):
    req = urllib.request.Request(url, headers={"Range": "bytes=0-0", "User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=120) as r:
        cr = r.headers.get("Content-Range", "")
        return int(cr.split("/")[1]) if "/" in cr else None

def fetch_blob(url, fbeg, fend_excl, total=None):
    if total is not None and fend_excl > total:
        fend_excl = total
    req = urllib.request.Request(url, headers={"Range": f"bytes={fbeg}-{fend_excl-1}", "User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=300) as r:
        if r.status == 200 and fbeg != 0:
            raise RuntimeError("server ignored Range")
        return r.read()

def fetch_region_rows(url, tbi, chrom, beg, end):
    found, ranges = region_file_ranges(tbi, chrom, beg, end)
    if not found:
        return []
    total = file_size(url)
    parts = []
    for fbeg, fend in ranges:
        raw = fetch_blob(url, fbeg, fend + 65536, total=total)
        parts.append(decompress_bgzf(raw))
    return b"".join(parts).split(b"\n")

def cache_one(ds_label, qtd, qts, chrom, beg, end):
    out = os.path.join(CACHE_DIR, f"{qtd}.{chrom}.tsv")
    if os.path.exists(out) and os.path.getsize(out) > 1000:
        print(f"[skip] {out} exists", flush=True)
        return out, None
    tbi = parse_tbi(tbi_path(qtd, qts))
    url = f"{BASE}/{qts}/{qtd}/{qtd}.all.tsv.gz"
    t0 = time.time()
    rows = fetch_region_rows(url, tbi, chrom, beg, end)
    with open(out, "wb") as f:
        for r in rows:
            if r.strip():
                f.write(r + b"\n")
    print(f"[OK] {qtd} {chrom}:{beg}-{end} rows={len(rows)} {time.time()-t0:.0f}s -> {out}", flush=True)
    return out, len(rows)

def job(ds):
    label, qtd, qts = ds
    for chrom, beg, end in REGIONS:
        cache_one(label, qtd, qts, chrom, beg, end)

if __name__ == "__main__":
    t_start = time.time()
    with cf.ThreadPoolExecutor(max_workers=5) as ex:
        futs = [ex.submit(job, d) for d in DATASETS]
        for f in cf.as_completed(futs):
            f.result()  # 抛出异常即失败
    print(f"ALL DONE in {time.time()-t_start:.0f}s; cache dir: {CACHE_DIR}")
    for fn in sorted(os.listdir(CACHE_DIR)):
        print("  ", fn, os.path.getsize(os.path.join(CACHE_DIR, fn)))
