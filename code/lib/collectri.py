# -*- coding: utf-8 -*-
"""code/lib/collectri.py —— CollecTRI（OmniPath）TF–靶基因网络的获取与版本校验

用途
----
`code/sc_regulatory/step0_build_cache.py` 与 `code/sc_regulatory/step2_regulon.py`
都需要 CollecTRI 的 TF→靶基因网络。该文件属**第三方原始资源**，按 `README.md`
「License」段的声明**不随本仓库重分发** ⇒ 第三方克隆本仓库后该文件必然缺失。
本模块提供 fallback 自动获取，使两脚本开箱即用（"文件存在时行为与原先完全一致"）。

获取顺序（取第一个命中者）
--------------------------
  1. 环境变量 `$IGAN_COLLECTRI` 指向的文件（显式指定，优先级最高）
  2. `<REPO>/data/sc_regulatory/CollecTRI_omnipath.tsv`（分析树内既有位置）
  3. `<REPO>/data/sc_regulatory/raw_cache/CollecTRI_omnipath.tsv`（本模块下载落点）
  4. `$RAWDATA/CollecTRI/collectri_omnipath.tsv`（手工下载到外部原始数据根）
  5. 从 OmniPath 自动下载 → 写入 (3)

⚠ 版本告警（重要，勿静默）
--------------------------
CollecTRI 是**持续更新的活资源**，OmniPath 不提供"按版本固定下载"的入口，
因此**在线查询结果会随时间变化**。本仓库分析所用的网络为
**完整的 CollecTRI 网络：64,516 条记录**（规范化 md5
`9296d2eb6f687dafe4067c233d5330d4`；1,201 个调控子、6,628 个不同靶基因），
它在 `step2_regulon.py` 的 `min_targets=5 / max_targets=None`（**不设靶数上限**）
口径下、**限定在表达面板（8,248 基因）内**产出
**778 个合格 regulon**，是论文 Table S10 / Fig. S6 的输入。

> 历史沿革：2026-10-01 之前，本包所用的是同一输入的**中断下载副本**
> （16,236 条完整记录，规范化 md5 `857151b015c24a1fe180bcc323162e46`），
> 在 `max_targets=500` 口径下只得 **317 个 regulon**。该副本已被完整网络取代，
> 去标识见 `README.md`「CollecTRI 版本说明」一节。

⇒ **用在线版本重跑，未必得到与论文冻结结果（Table S10 / Fig. S6）完全一致的数值**
（网络自 2026-10-01 之后仍会继续更新）。
本模块在记录数/校验和与该完整网络不符时**打印显著告警**（但不中断），由使用者判断；
告警不会静默通过。

数据格式
--------
OmniPath `interactions` 接口 TSV，关键列为
`source` / `target`（UniProt）与 `source_genesymbol` / `target_genesymbol`。
注意 `genesymbols=1` 是必需的：不带该参数时接口**不返回** `*_genesymbol` 两列。
"""
import gzip
import hashlib
import os
import sys
import urllib.request

# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))          # <REPO>/code/lib
_CODE = os.path.dirname(_HERE)                              # <REPO>/code
if _CODE not in sys.path:
    sys.path.insert(0, _CODE)
import _paths  # noqa: E402

# --------------------------------------------------------------------------
COLLECTRI_URL = ("https://omnipathdb.org/interactions"
                 "?datasets=collectri&genesymbols=1&format=tsv")
FILENAME = "CollecTRI_omnipath.tsv"
REQUIRED_COLS = ("source_genesymbol", "target_genesymbol")

# 本仓库分析所用网络的指纹（规范化字节：表头 + 全部完整记录，末尾补 1 个换行）
SNAPSHOT_RECORDS = 64516
SNAPSHOT_MD5 = "9296d2eb6f687dafe4067c233d5330d4"
SNAPSHOT_REGULONS = 778          # 表达面板内 / min_targets=5 / max_targets=None（无上限）下的合格 regulon 数
_UA = {"User-Agent": "Mozilla/5.0 (research-repro; CollecTRI fetch)",
       "Accept-Encoding": "gzip"}


def _split_records(text):
    """切出「表头 + 完整数据行」。末尾若不完整（字段数 < 表头列数）则丢弃。"""
    lines = text.split("\n")
    while lines and lines[-1] == "":
        lines.pop()
    if not lines:
        return []
    ncol = len(lines[0].split("\t"))
    if len(lines[-1].split("\t")) != ncol:
        lines = lines[:-1]
    return lines


def summarise(path):
    """→ dict(n_records, columns, md5_norm)；md5_norm 对「表头+完整记录+末尾换行」计算。"""
    raw = open(path, "rb").read()
    lines = _split_records(raw.decode("utf-8", "replace"))
    if not lines:
        return {"n_records": 0, "columns": [], "md5_norm": hashlib.md5(b"").hexdigest()}
    norm = ("\n".join(lines) + "\n").encode("utf-8")
    return {"n_records": len(lines) - 1,
            "columns": lines[0].split("\t"),
            "md5_norm": hashlib.md5(norm).hexdigest()}


def _download(dst, tries=3):
    """从 OmniPath 下载并落盘；返回落盘字节数。失败重试，最终抛出并给出可读原因。"""
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    last = None
    for k in range(1, tries + 1):
        try:
            print(f"[collectri] 下载 {COLLECTRI_URL}  (第 {k}/{tries} 次)", flush=True)
            req = urllib.request.Request(COLLECTRI_URL, headers=_UA)
            with urllib.request.urlopen(req, timeout=600) as r:
                raw = r.read()
            if raw[:2] == b"\x1f\x8b":                 # 服务器按请求返回 gzip
                raw = gzip.decompress(raw)
            text = raw.decode("utf-8")
            lines = _split_records(text)
            if len(lines) < 2 or not all(c in lines[0] for c in REQUIRED_COLS):
                raise RuntimeError(
                    "返回内容不含所需列 " + "/".join(REQUIRED_COLS)
                    + f"（实际表头：{lines[0][:120] if lines else '空'}）")
            tmp = dst + ".part"
            with open(tmp, "w", encoding="utf-8", newline="") as f:
                f.write("\n".join(lines) + "\n")
            os.replace(tmp, dst)
            return os.path.getsize(dst)
        except Exception as e:                          # noqa: BLE001
            last = e
            print(f"[collectri] 第 {k} 次失败：{type(e).__name__}: {e}", flush=True)
    raise RuntimeError(
        f"无法从 OmniPath 获取 CollecTRI：{last}\n"
        f"  可手工下载后放到 {dst}\n"
        f"  或把已有文件路径写入环境变量 IGAN_COLLECTRI 后重跑。")


def warn_if_not_snapshot(info, path):
    """把版本差异显著打印出来（不中断）。"""
    n, h = info["n_records"], info["md5_norm"]
    if h == SNAPSHOT_MD5:
        print(f"[collectri] 版本校验：与论文所用完整网络一致（{n:,} 条记录）", flush=True)
        return True
    print("=" * 78, flush=True)
    print("[collectri] ⚠ 版本告警：当前 CollecTRI 与论文所用网络不同", flush=True)
    print(f"  文件            : {path}", flush=True)
    print(f"  论文所用网络    : {SNAPSHOT_RECORDS:,} 条记录 · md5 {SNAPSHOT_MD5}"
          f"  →  合格 regulon {SNAPSHOT_REGULONS} 个（min_targets=5，无上限）", flush=True)
    print(f"  当前获取        : {n:,} 条记录 · md5 {h}", flush=True)
    print("  影响            : min_targets=5 / max_targets=None 口径下 regulon 集合会变化", flush=True)
    print("  ⇒ 用当前版本重跑，Table S10 / Fig. S6 的数值不会与论文冻结结果完全一致。", flush=True)
    print("    要严格复现论文数值：使用论文所用完整网络（64,516 条完整记录）。", flush=True)
    print("=" * 78, flush=True)
    return False


def resolve(path=None, data_dir=None, rawdata_dir=None, download=True, verbose=True):
    """返回可用的 CollecTRI 文件路径（必要时下载）。"""
    if path:
        if not os.path.isfile(path):
            raise FileNotFoundError(f"[collectri] 指定文件不存在：{path}")
        if verbose:
            warn_if_not_snapshot(summarise(path), path)
        return path

    env = os.environ.get("IGAN_COLLECTRI")
    if env:
        if not os.path.isfile(env):
            raise FileNotFoundError(f"[collectri] $IGAN_COLLECTRI 指向的文件不存在：{env}")
        if verbose:
            warn_if_not_snapshot(summarise(env), env)
        return env

    data_dir = data_dir or _paths.data("sc_regulatory")
    rawdata_dir = rawdata_dir or _paths.RAWDATA
    cache = os.path.join(data_dir, "raw_cache", FILENAME)
    cands = [
        os.path.join(data_dir, FILENAME),
        cache,
        os.path.join(rawdata_dir, "CollecTRI", "collectri_omnipath.tsv"),
    ]
    for c in cands:
        if os.path.isfile(c):
            if verbose:
                warn_if_not_snapshot(summarise(c), c)
            return c

    if not download:
        raise FileNotFoundError(
            "[collectri] 未找到 CollecTRI 文件，且 download=False。查找过：\n  "
            + "\n  ".join(cands))

    size = _download(cache)
    print(f"[collectri] 已落盘 {cache}（{size/1024/1024:.1f} MB）", flush=True)
    if verbose:
        warn_if_not_snapshot(summarise(cache), cache)
    return cache


def load(path=None, **kw):
    """→ pandas.DataFrame：完整 CollecTRI 表（供 step0/step2 直接使用）。"""
    import pandas as pd
    p = resolve(path, **kw)
    df = pd.read_csv(p, sep="\t", dtype=str)
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise RuntimeError(f"[collectri] {p} 缺列：{missing}（实际列：{list(df.columns)}）")
    return df


if __name__ == "__main__":
    p = resolve()
    info = summarise(p)
    print(f"path       : {p}")
    print(f"records    : {info['n_records']:,}")
    print(f"md5_norm   : {info['md5_norm']}")
    print(f"snapshot?  : {info['md5_norm'] == SNAPSHOT_MD5}")
    print(f"columns    : {info['columns']}")
