# -*- coding: utf-8 -*-
"""
_tmp_patch_o_drop_stray_allprobe.py   (2026-09-10, 作者裁定「删除」)

作者裁定：删除 submission/additional_file_1/Figure_S1_GSE73953_IgAN_vs_HC_allprobe.csv
—— 该件为逐探针差异明细（888 B，表头 + 7 探针行，覆盖 5 个轴基因），前缀 `Figure_S1_`
易被误认为图件；未被任何交付图注 / 补充材料索引 / 渲染脚本引用。

删除依据（删前核验）：
  1. md5 三副本全同 = 789314d5c6b21d33e2a561070c105dda
       - 权威源（保留不删）：阶段1/M1_GEO表达锚点/GSE73953_IgAN_vs_HC_allprobe.csv
       - 复现镜像（随删）：阶段4.5_复现包/reproducibility/code/additional_file_1/Figure_S1_GSE73953_IgAN_vs_HC_allprobe.csv
       - 交付件（本例删除）：阶段4.5_复现包/submission/additional_file_1/Figure_S1_GSE73953_IgAN_vs_HC_allprobe.csv
  2. S1 渲染脚本 figure_s1_3datasets_boxplot.py 只读 _groups.tsv / _5genes_expr.csv / _best.csv，
     **不读** _allprobe.csv —— 删除不影响 S1 复现。
  3. 交付图注 S1 的 per-gene logFC / P 值指向 **Table 2**，未引用该 CSV。

恢复方式（如需）：cp 阶段1/M1_GEO表达锚点/GSE73953_IgAN_vs_HC_allprobe.csv <目标>
不涉及任何数值、结论、图件或图注改动。
"""
# --- 统一路径解析（开箱即用）：定位 code/_paths.py -------------------------
import os as _os, sys as _sys
_C = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _C if _os.path.isfile(_os.path.join(_C, "_paths.py"))
                else _os.path.dirname(_C))
import _paths  # noqa: E402  （复现包路径解析；IGAN_* 环境变量可覆盖根目录）
# --------------------------------------------------------------------------
import os
import hashlib

ROOT = _paths.LEGACY

TARGETS = [
    r"阶段4.5_复现包\submission\additional_file_1\Figure_S1_GSE73953_IgAN_vs_HC_allprobe.csv",
    r"阶段4.5_复现包\reproducibility\code\additional_file_1\Figure_S1_GSE73953_IgAN_vs_HC_allprobe.csv",
]

CANONICAL = r"阶段1\M1_GEO表达锚点\GSE73953_IgAN_vs_HC_allprobe.csv"
EXPECT_MD5 = "789314d5c6b21d33e2a561070c105dda"


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    # 0) 权威源必须存在且与预期 md5 一致
    cp = ROOT(CANONICAL)
    if not os.path.exists(cp):
        raise SystemExit("ABORT: 权威源缺失，拒绝删除 -> " + cp)
    cm = md5(cp)
    if cm != EXPECT_MD5:
        raise SystemExit("ABORT: 权威源 md5 异常 %s != %s" % (cm, EXPECT_MD5))
    print("[keep] 权威源 %s\n        md5=%s" % (CANONICAL, cm))

    removed = 0
    for rel in TARGETS:
        p = ROOT(rel)
        if not os.path.exists(p):
            print("[skip] 已不存在: %s" % rel)
            continue
        m = md5(p)
        if m != EXPECT_MD5:
            print("[WARN] md5 不符，跳过: %s (got %s)" % (rel, m))
            continue
        os.remove(p)
        removed += 1
        print("[del ] %s (md5=%s, 与权威源一致)" % (rel, m))

    print("\n完成：删除 %d 件；权威源保留，随时可恢复。" % removed)


if __name__ == "__main__":
    main()
