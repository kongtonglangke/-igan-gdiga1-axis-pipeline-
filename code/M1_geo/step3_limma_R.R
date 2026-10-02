# -*- coding: utf-8 -*-
# ============================================================
#  M1 · 步骤 1.1 GEO 表达锚点分析 (R+limma 等价脚本)
#  ----------------------------------------------------------
#  本脚本与 step1_parse_and_de.py (Python 沙箱执行版) 在设计上等价:
#    1) 同样解析 3 个 GEO series matrix (GSE73953 / GSE115857 / GSE93798);
#    2) 同样用平台注释文件把探针映射到 5 个机制轴基因 (C1GALT1, C1GALT1C1,
#       ST6GALNAC2, GALNT2, GALNT12);
#    3) 差异分析采用 limma::lmFit + eBayes + topTable (与 Python Welch t-test
#       语义近似; 沙箱执行以 Python 实现, 本脚本供 R 环境复现 / 报告入档);
#    4) 输出差异表 (best.csv / allprobe.csv)、分组 (groups.tsv)、
#       5 基因表达子集 (5genes_expr.csv, log2 尺度) 与 ggplot2 箱线图。
#
#  平台注释依赖 (请在 R 环境中预先安装):
#    - GPL4133: 直接读 Agilent 注释 (从 .annot.gz 解压后读取)
#    - GPL14951: Bioconductor `illuminaHumanv4.db_1.26.0` (ILMN -> Entrez)
#    - GPL22945: BrainArray custom CDF v19 (Entrez-on-probe) -- 不需额外包,
#                探针即为 {EntrezGeneID}_at
#
#  备注:
#   * 沙箱无 R, 实际差异分析已由 step1_parse_and_de.py 跑通; 输出一致
#     参见 {GSE}_*_best.csv 与 figures/boxplot_*.png.
#   * 本脚本与 Python 版共用数据路径约定 (RAW=/00_rawdata/GEO_M1,
#     OUT=./M1_GEO表达锚点), 在 Windows / Linux / macOS 上 R 4.x + Bioc 3.16+
#     可直接 source() 执行。
# ============================================================
suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  for (p in c("GEOquery","limma","ggplot2","dplyr","tidyr","readr",
              "illuminaHumanv4.db")) {
    if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE)
  }
  library(GEOquery); library(limma); library(ggplot2)
  library(dplyr); library(tidyr); library(readr)
  library(illuminaHumanv4.db)
})

# ---------- 路径约定 (与 Python 版一致) ----------
# 由复现包 code/_paths.py 同源解析；可用 IGAN_REPRO_ROOT / IGAN_RAWDATA_ROOT 覆盖
.r_args   <- commandArgs(trailingOnly = FALSE)
.r_script <- sub("^--file=", "", .r_args[grep("^--file=", .r_args)])
.repro    <- Sys.getenv("IGAN_REPRO_ROOT",
                        unset = normalizePath(file.path(dirname(.r_script), ".."),
                                              mustWork = FALSE))
.rawdata  <- Sys.getenv("IGAN_RAWDATA_ROOT",
                        unset = normalizePath(file.path(.repro, "..", "00_rawdata"),
                                              mustWork = FALSE))
RAW <- normalizePath(file.path(.rawdata, "GEO_M1"), mustWork = FALSE)
OUT <- file.path(.repro, "data", "M1_geo")
dir.create(file.path(OUT, "figures"), showWarnings = FALSE, recursive = TRUE)

# ---------- 5 机制轴基因 (symbol -> Entrez) ----------
AXIS <- c("C1GALT1"=56913, "C1GALT1C1"=29071, "ST6GALNAC2"=10610,
          "GALNT2"=2590, "GALNT12"=79695)

# ---------- 系列矩阵解析 + 分组 (与 Python 版一致) ----------
parse_series_matrix <- function(gz_path) {
  # 简化: 用 GEOquery::getGEO 拉软格式; 离线模式可用本地 series matrix
  if (Sys.which("curl") != "" && RCurl::url.exists(paste0("ftp://ftp.ncbi.nlm.nih.gov/geo/series/",
                                                          substr(basename(gz_path), 1, nchar(basename(gz_path))-3)))) {
    # 优先尝试 GEOquery 在线 (用户 R 联网时)
    gse_id <- sub("_series_matrix.*", "", basename(gz_path))
    gse <- tryCatch(getGEO(gse_id, destdir = file.path(RAW, "geo_cache")),
                    error = function(e) NULL)
    if (!is.null(gse)) {
      eset <- gse[[1]]
      pheno <- pData(eset); mat <- exprs(eset)
      return(list(expr = mat, pheno = pheno))
    }
  }
  # 离线 fallback: 解本地 series matrix
  con <- gzfile(gz_path, "rt"); lines <- readLines(con, warn = FALSE); close(con)
  body <- lines[grep("!series_matrix_table_begin", lines):
                grep("!series_matrix_table_end", lines)]
  cols <- strsplit(sub("^\"", "", sub("\"$", "", body[1])), "\t")[[1]]
  df <- read.table(text = paste(body[-1], collapse = "\n"),
                   sep = "\t", header = FALSE, fill = TRUE, quote = "\"")
  colnames(df) <- cols
  rownames(df) <- df[[1]]; df[[1]] <- NULL
  df <- as.matrix(data.frame(lapply(df, as.numeric), row.names = rownames(df),
                             check.names = FALSE))
  pheno <- data.frame(
    sample_id = cols[-1],
    title = sapply(strsplit(lines[grep("!Sample_title", lines)], "\t"), `[`, 2),
    source = sapply(strsplit(lines[grep("!Sample_source_name_ch1", lines)], "\t"), `[`, 2),
    stringsAsFactors = FALSE
  )
  list(expr = df, pheno = pheno)
}

# ---------- 分组规则 (与 Python 一致) ----------
make_groups <- function(pheno, gse) {
  acc <- pheno$sample_id
  if (gse == "GSE73953") {
    title <- tolower(pheno$title)
    grp <- ifelse(grepl("iga nephrop|igan", title), "IgAN",
           ifelse(grepl("membranous", title), "MN",
           ifelse(grepl("healthy|control", title), "HC", "OTHER")))
  } else if (gse == "GSE115857") {
    src <- tolower(pheno$source)
    grp <- ifelse(grepl("igan", src), "IgAN",
           ifelse(grepl("living donor", src), "Ctrl_LD",
           ifelse(grepl("membranous", src), "MN",
           ifelse(grepl("minimal change", src), "MCD",
           ifelse(grepl("focal", src), "FSGS", "OTHER")))))
  } else { # GSE93798
    src <- tolower(pheno$source)
    grp <- ifelse(grepl("control|^con$", src), "Ctrl",
           ifelse(grepl("igan", src), "IgAN", "OTHER"))
  }
  setNames(grp, acc)
}

# ---------- 平台 → 探针→基因映射 (与 Python 一致) ----------
map_probes <- function(gse, expr_mat) {
  probes <- rownames(expr_mat)
  if (gse == "GSE93798") {
    # 探针 = {EntrezGeneID}_at, 直接匹配 AXIS
    m <- setNames(rep(NA_character_, length(probes)), probes)
    for (sym in names(AXIS)) {
      m[probes == paste0(AXIS[sym], "_at")] <- sym
    }
    return(m[!is.na(m)])
  }
  if (gse == "GSE73953") {
    # GPL4133.annot: ID Gene symbol 两列
    ann <- read.delim(gzfile(file.path(RAW, "GPL4133.annot.gz")),
                      comment.char = "!", stringsAsFactors = FALSE)
    m <- setNames(ann$Gene.symbol, as.character(ann$ID))
    m <- m[m %in% names(AXIS)]
    return(m)
  }
  if (gse == "GSE115857") {
    # illuminaHumanv4.db: ILMN_ -> Entrez (通过 probes() 选择列)
    ilmn <- toTable(illuminaHumanv4PROBE)  # probe_id, gene_id
    keep <- ilmn$gene_id %in% as.character(AXIS)
    m <- setNames(names(AXIS)[match(ilmn$gene_id[keep], as.character(AXIS))],
                  ilmn$probe_id[keep])
    return(m)
  }
  character(0)
}

# ---------- limma 差异分析 ----------
run_limma <- function(gse, gz_path, case, ctrl, log2 = FALSE) {
  cat(sprintf("\n== %s (%s vs %s) ==\n", gse, case, ctrl))
  parsed <- parse_series_matrix(gz_path)
  expr <- parsed$expr
  if (log2) expr <- log2(pmax(expr, 1e-3))
  pheno <- parsed$pheno
  groups <- make_groups(pheno, gse)
  sample_keep <- names(groups)[groups %in% c(case, ctrl)]
  expr <- expr[, sample_keep, drop = FALSE]
  groups <- groups[sample_keep]
  # 平台映射
  probe_gene <- map_probes(gse, expr)
  if (length(probe_gene) == 0) { cat("  WARN no probe mapped\n"); return(NULL) }
  # limma
  design <- model.matrix(~ 0 + factor(groups, levels = c(ctrl, case)))
  colnames(design) <- c(ctrl, case)
  fit <- lmFit(expr[rownames(expr) %in% names(probe_gene), , drop = FALSE], design)
  contrast <- makeContrasts(IgAN_vs_Ctrl = paste0(case, "-", ctrl), levels = design)
  fit2 <- eBayes(contrasts.fit(fit, contrast))
  # 全部探针结果
  allp <- topTable(fit2, coef = 1, number = Inf, sort.by = "none")
  allp$probe <- rownames(allp); allp$gene <- probe_gene[allp$probe]
  allp$n_case <- sum(groups == case); allp$n_ctrl <- sum(groups == ctrl)
  # 多探针取最显著
  best <- allp %>% group_by(gene) %>% slice_min(P.Value, n = 1) %>% ungroup()
  cat(sprintf("  probe mapped: %d | genes: %s\n", nrow(allp), paste(sort(unique(best$gene)), collapse = ", ")))
  print(best[, c("gene","probe","n_case","n_ctrl","logFC","t","P.Value")])
  # 输出
  tag <- paste0(case, "_vs_", ctrl)
  write.csv(best, file.path(OUT, sprintf("%s_%s_best.csv", gse, tag)), row.names = FALSE)
  write.csv(allp, file.path(OUT, sprintf("%s_%s_allprobe.csv", gse, tag)), row.names = FALSE)
  write.table(data.frame(sample = names(groups), group = groups, stringsAsFactors = FALSE),
              file.path(OUT, sprintf("%s_groups.tsv", gse)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  # 5 基因表达子集 (多探针取均值)
  keep <- allp$probe
  sub <- expr[keep, , drop = FALSE]
  gmap <- setNames(allp$gene, allp$probe)
  rownames(sub) <- gmap[rownames(sub)]
  agg <- rowsum(t(sub), group = rownames(sub)) / as.numeric(table(rownames(sub)))
  agg <- t(agg)  # gene x sample
  write.csv(agg, file.path(OUT, sprintf("%s_5genes_expr.csv", gse)))
  invisible(list(expr = agg, groups = groups, best = best))
}

# ---------- ggplot2 箱线图 ----------
plot_boxplot <- function(gse, title, case, ctrl, agg, groups, best) {
  df <- as.data.frame(agg) %>% mutate(gene = rownames(agg)) %>%
    pivot_longer(-gene, names_to = "sample", values_to = "expr")
  df$group <- groups[df$sample]
  df <- merge(df, best[, c("gene","logFC","P.Value")], by = "gene", all.x = TRUE)
  df$gene <- factor(df$gene, levels = names(AXIS))
  p <- ggplot(df, aes(group, expr, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.5) +
    geom_jitter(width = 0.06, size = 0.7, alpha = 0.45) +
    scale_fill_manual(values = c("IgAN" = "#d62728", "HC" = "#1f77b4",
                                 "Ctrl" = "#1f77b4", "Ctrl_LD" = "#1f77b4")) +
    facet_wrap(~ gene, scales = "free_y", nrow = 1) +
    labs(title = sprintf("%s — %s", gse, title), x = NULL, y = "expression (log2)") +
    theme_bw(base_family = "sans") + theme(legend.position = "none")
  ggsave(file.path(OUT, "figures", sprintf("boxplot_%s.png", gse)),
         p, width = 14, height = 3.2, dpi = 200)
  ggsave(file.path(OUT, "figures", sprintf("boxplot_%s.pdf", gse)),
         p, width = 14, height = 3.2)
}

# ---------- 主流程 ----------
main <- function() {
  r1 <- run_limma("GSE73953", file.path(RAW, "GSE73953_series_matrix.txt.gz"),
                  "IgAN", "HC", log2 = TRUE)
  r2 <- run_limma("GSE115857", file.path(RAW, "GSE115857_series_matrix.txt.gz"),
                  "IgAN", "Ctrl_LD")
  r3 <- run_limma("GSE93798", file.path(RAW, "GSE93798_series_matrix.txt.gz"),
                  "IgAN", "Ctrl")
  if (!is.null(r1)) plot_boxplot("GSE73953", "PBMC · IgAN vs HC-pooled", "IgAN", "HC", r1$expr, r1$groups, r1$best)
  if (!is.null(r2)) plot_boxplot("GSE115857", "Kidney biopsy bulk · IgAN vs Living donor", "IgAN", "Ctrl_LD", r2$expr, r2$groups, r2$best)
  if (!is.null(r3)) plot_boxplot("GSE93798", "Glomeruli · IgAN vs Control", "IgAN", "Ctrl", r3$expr, r3$groups, r3$best)
}

if (sys.nframe() == 0) main()
