## =====================================================================
## 03c_coloc_mqtl_eqtl.R  -  mQTL x C1GALT1 cis-eQTL colocalisation
##
## Replaces the pure-Python ABF used for Table S5d with coloc::coloc.abf, so that
## every region-level colocalisation in the paper is computed by the official
## coloc package.
##
## in : data/mqtl/cpg_full_assoc_<CpG>.tsv   (GoDMC; fixed-effects `se`, NOT `se_mre`)
##      data/eqtl_for_coloc/QTD000xxx.7.tsv  (eQTL Catalogue r8, C1GALT1 on chr7)
## out: out/14_coloc_mqtl_eqtl_R.tsv
## =====================================================================
.this_dir <- dirname(normalizePath(sub("^--file=", "",
  commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]),
  winslash = "/", mustWork = FALSE))
source(file.path(.this_dir, "00_config.R"))

suppressPackageStartupMessages(library(coloc))

ENSG_C1GALT1 <- "ENSG00000106392"
QTDS <- c(QTD000356 = "GTEx_v10_blood", QTD000261 = "GTEx_v10_kidney_cortex",
          QTD000608 = "OneK1K_B_naive", QTD000607 = "OneK1K_B_memory",
          QTD000606 = "OneK1K_B_intermediate")

read_mqtl <- function(cpg) {
  f <- file.path(MQLT_D, sprintf("cpg_full_assoc_%s.tsv", cpg))
  d <- data.table::fread(f, sep = "\t", header = TRUE, showProgress = FALSE)
  d <- d[, .(rsid, beta = as.numeric(beta_are_a1), se = as.numeric(se),
             a1 = toupper(a1), a2 = toupper(a2))]
  d <- d[!is.na(rsid) & rsid != "" & !is.na(beta) & !is.na(se) & se > 0]
  d[, .SD[!duplicated(rsid)]]
}

read_eqtl7 <- function(qtd) {
  f <- file.path(ROOT, "data/eqtl_for_coloc", sprintf("%s.7.tsv", qtd))
  if (!file.exists(f)) return(NULL)
  d <- data.table::fread(f, sep = "\t", header = TRUE, showProgress = FALSE)
  d <- d[gene == ENSG_C1GALT1]
  if (nrow(d) == 0) return(NULL)
  d[, `:=`(rsid = rsid, alt = toupper(alt), ref = toupper(ref))]
  d <- d[!is.na(rsid) & rsid != "" & !is.na(beta) & !is.na(se) & se > 0]
  d[, .SD[!duplicated(rsid)]]
}

rows <- list()
for (cpg in CPG_INFO$cpg) {
  mq <- read_mqtl(cpg)
  for (q in names(QTDS)) {
    eq <- read_eqtl7(q)
    if (is.null(eq)) {
      rows[[length(rows) + 1]] <- data.frame(
        cpg = cpg, qtd = q, context = QTDS[[q]], n_shared = 0L,
        mqtl_top = NA_character_, eqtl_top = NA_character_, eqtl_top_P = NA_real_,
        PPH0 = NA_real_, PPH1 = NA_real_, PPH2 = NA_real_, PPH3 = NA_real_, PPH4 = NA_real_,
        note = "not computable: no cis-eQTL record")
      next
    }
    m <- merge(mq, eq, by = "rsid", suffixes = c(".m", ".e"))
    ## mQTL effect allele = a1; eQTL effect allele = alt
    same <- m$a1 == m$alt & m$a2 == m$ref
    flip <- m$a1 == m$ref & m$a2 == m$alt
    keep <- same | flip
    m <- m[keep, , drop = FALSE]
    if (nrow(m) < 50) {
      rows[[length(rows) + 1]] <- data.frame(
        cpg = cpg, qtd = q, context = QTDS[[q]], n_shared = nrow(m),
        mqtl_top = NA_character_, eqtl_top = NA_character_, eqtl_top_P = NA_real_,
        PPH0 = NA_real_, PPH1 = NA_real_, PPH2 = NA_real_, PPH3 = NA_real_, PPH4 = NA_real_,
        note = sprintf("not computable: only %d shared variants", nrow(m)))
      next
    }
    b <- as.numeric(m[["beta.e"]])
    f <- as.logical(flip[keep])
    b[f] <- -b[f]
    m[["beta.e"]] <- b

    d1 <- list(snp = m$rsid, beta = as.numeric(m[["beta.m"]]), varbeta = as.numeric(m[["se.m"]])^2,
               position = m$pos, type = "quant", N = 25000, sdY = 1)
    d2 <- list(snp = m$rsid, beta = m[["beta.e"]], varbeta = as.numeric(m[["se.e"]])^2,
               position = m$pos, type = "quant", N = 700, sdY = 1)
    r <- try(coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5), silent = TRUE)
    if (inherits(r, "try-error")) {
      rows[[length(rows) + 1]] <- data.frame(
        cpg = cpg, qtd = q, context = QTDS[[q]], n_shared = nrow(m),
        mqtl_top = NA_character_, eqtl_top = NA_character_, eqtl_top_P = NA_real_,
        PPH0 = NA_real_, PPH1 = NA_real_, PPH2 = NA_real_, PPH3 = NA_real_, PPH4 = NA_real_,
        note = "coloc.abf failed")
      next
    }
    s <- r$summary
    rows[[length(rows) + 1]] <- data.frame(
      cpg = cpg, qtd = q, context = QTDS[[q]], n_shared = as.integer(s[["nsnps"]]),
      mqtl_top = m$rsid[which.max(abs(m[["beta.m"]] / m[["se.m"]]))],
      eqtl_top = m$rsid[which.min(eq$se[match(m$rsid, eq$rsid)]^2)],
      eqtl_top_P = NA_real_,
      PPH0 = s[["PP.H0.abf"]], PPH1 = s[["PP.H1.abf"]], PPH2 = s[["PP.H2.abf"]],
      PPH3 = s[["PP.H3.abf"]], PPH4 = s[["PP.H4.abf"]],
      note = sprintf("coloc.abf; %d shared variants, %d flipped", sum(f), sum(f)))
  }
}

out <- do.call(rbind, rows)
data.table::fwrite(out, file.path(DIR_OUT, "14_coloc_mqtl_eqtl_R.tsv"), sep = "\t")
cat("[03c] comparisons:", nrow(out), " computable:", sum(!is.na(out$PPH4)), "\n")
cat("[03c] max PPH4 =", round(max(out$PPH4, na.rm = TRUE), 4), "\n")
