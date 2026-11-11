## =====================================================================
## 03b_coloc_eqtl_igan.R  -  eQTL x IgAN colocalisation, 5 genes x 5 contexts
##
## Replaces the pure-Python ABF implementation used up to 2026-11-11 with the
## official coloc package, so that every region-level colocalisation in the paper
## is computed by coloc::coloc.abf.
##
## in : data/igan_regions_for_coloc/chr{1,7,9,17}.tsv   (IgAN extracts, see
##      _extract_igan_regions.py; chrX is absent from GCST90018866)
##      $EQTL_CACHE/QTD000xxx.<chr>.tsv                  (eQTL Catalogue r8, long format)
## out: out/13_coloc_eqtl_igan_R.tsv
##
## Priors are the pre-specified ones: p1 = p2 = 1e-4, p12 = 1e-5.
## The IgAN dataset is analysed as case-control with its own case fraction; the
## eQTL side is quantitative. sd.prior is left at the coloc default for both types.
## =====================================================================
.this_dir <- dirname(normalizePath(sub("^--file=", "",
  commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]),
  winslash = "/", mustWork = FALSE))
source(file.path(.this_dir, "00_config.R"))

suppressPackageStartupMessages(library(coloc))

IGAN_D <- file.path(ROOT, "data/igan_regions_for_coloc")

AXIS <- data.frame(
  gene = c("C1GALT1", "C1GALT1C1", "GALNT12", "GALNT2", "ST6GALNAC2"),
  ensg = c("ENSG00000106392", "ENSG00000171155", "ENSG00000119514",
           "ENSG00000143641", "ENSG00000070731"),
  chr  = c("7", "X", "9", "1", "17"),
  start = c(7156934, 120614846, 98807613, 230057990, 76564424),
  end   = c(7259333, 120630115, 98851267, 230282128, 76586956),
  stringsAsFactors = FALSE)

QTDS <- c(QTD000356 = "GTEx_v10_blood", QTD000261 = "GTEx_v10_kidney_cortex",
          QTD000608 = "OneK1K_B_naive", QTD000607 = "OneK1K_B_memory",
          QTD000606 = "OneK1K_B_intermediate")

FLANK <- 1e6
S_IGAN <- N_IGAN_CASE / N_IGAN

## ---- eQTL Catalogue long format: no header, columns taken by position -----
## 1 gene_id | 2 chr | 3 pos | 4 ref | 5 alt | 6 variant | 8 maf | 9 p | 10 beta | 11 se | 19 rsid
read_eqtl <- function(qtd, chr_, ensg, lo, hi) {
  f <- file.path(ROOT, "data/eqtl_for_coloc", sprintf("%s.%s.tsv", qtd, chr_))
  if (!file.exists(f)) return(NULL)
  d <- data.table::fread(f, sep = "\t", header = TRUE, showProgress = FALSE)
  if (nrow(d) == 0) return(NULL)
  d <- d[gene == ensg & !is.na(rsid) & rsid != "" & pos >= lo & pos <= hi]
  if (nrow(d) == 0) return(NULL)
  d <- d[!is.na(beta) & !is.na(se) & se > 0]
  d[, .SD[!duplicated(rsid)]]
}

read_igan <- function(chr_) {
  f <- file.path(IGAN_D, sprintf("chr%s.tsv", chr_))
  if (!file.exists(f)) return(NULL)
  d <- data.table::fread(f, sep = "\t", header = TRUE, showProgress = FALSE)
  d[, `:=`(pos = as.numeric(pos), beta = as.numeric(beta), se = as.numeric(se))]
  d <- d[!is.na(rsid) & rsid != "" & !is.na(beta) & !is.na(se) & se > 0]
  d[, .SD[!duplicated(rsid)]]
}

## ---- allele alignment: flip the eQTL effect when it is reported on the other allele
align <- function(e, g) {
  m <- merge(e, g, by = "rsid", suffixes = c(".e", ".g"))
  if (nrow(m) == 0) return(m)
  ## eQTL Catalogue reports beta for the ALT allele, so ALT is the effect allele
  same <- (toupper(m$alt) == toupper(m$ea) & toupper(m$ref) == toupper(m$oa))
  flip <- (toupper(m$alt) == toupper(m$oa) & toupper(m$ref) == toupper(m$ea))
  keep <- same | flip
  m <- m[keep, , drop = FALSE]
  m$flip <- flip[keep]
  b <- as.numeric(m[["beta.e"]])
  f <- as.logical(m$flip)
  b[f] <- -b[f]
  m[["beta.e"]] <- b
  m
}

rows <- list()
for (i in seq_len(nrow(AXIS))) {
  g <- AXIS[i, ]
  if (g$chr == "X") {
    for (q in names(QTDS)) {
      rows[[length(rows) + 1]] <- data.frame(
        gene = g$gene, qtd = q, context = QTDS[[q]], n_snp = 0L,
        PPH0 = NA_real_, PPH1 = NA_real_, PPH2 = NA_real_, PPH3 = NA_real_, PPH4 = NA_real_,
        note = "not computable: C1GALT1C1 is on chrX and GCST90018866 has no X")
    }
    next
  }
  igan <- read_igan(g$chr)
  for (q in names(QTDS)) {
    e <- read_eqtl(q, g$chr, g$ensg, g$start - FLANK, g$end + FLANK)
    if (is.null(e) || nrow(e) == 0 || is.null(igan)) {
      rows[[length(rows) + 1]] <- data.frame(
        gene = g$gene, qtd = q, context = QTDS[[q]], n_snp = 0L,
        PPH0 = NA_real_, PPH1 = NA_real_, PPH2 = NA_real_, PPH3 = NA_real_, PPH4 = NA_real_,
        note = "not computable: no cis-eQTL record for this gene in this context")
      next
    }
    m <- align(e, igan)
    if (nrow(m) < 50) {
      rows[[length(rows) + 1]] <- data.frame(
        gene = g$gene, qtd = q, context = QTDS[[q]], n_snp = nrow(m),
        PPH0 = NA_real_, PPH1 = NA_real_, PPH2 = NA_real_, PPH3 = NA_real_, PPH4 = NA_real_,
        note = sprintf("not computable: only %d shared variants", nrow(m)))
      next
    }
    ## dataset1 = expression (H1 = eQTL only), dataset2 = IgAN (H2 = disease only),
    ## so the hypothesis labels match the manuscript
    d1 <- list(snp = m$rsid, beta = m$beta.e, varbeta = m$se.e^2,
               position = m$pos.e, type = "quant", N = 700, sdY = 1)
    d2 <- list(snp = m$rsid, beta = m$beta.g, varbeta = m$se.g^2,
               position = m$pos.e, type = "cc", s = S_IGAN, N = N_IGAN)
    r <- try(coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5), silent = TRUE)
    if (inherits(r, "try-error")) {
      rows[[length(rows) + 1]] <- data.frame(
        gene = g$gene, qtd = q, context = QTDS[[q]], n_snp = nrow(m),
        PPH0 = NA_real_, PPH1 = NA_real_, PPH2 = NA_real_, PPH3 = NA_real_, PPH4 = NA_real_,
        note = "coloc.abf failed")
      next
    }
    s <- r$summary
    rows[[length(rows) + 1]] <- data.frame(
      gene = g$gene, qtd = q, context = QTDS[[q]], n_snp = as.integer(s[["nsnps"]]),
      PPH0 = s[["PP.H0.abf"]], PPH1 = s[["PP.H1.abf"]], PPH2 = s[["PP.H2.abf"]],
      PPH3 = s[["PP.H3.abf"]], PPH4 = s[["PP.H4.abf"]],
      note = sprintf("coloc.abf; %d shared variants, %d allele-flipped",
                     sum(m$flip), sum(m$flip)))
  }
}

out <- do.call(rbind, rows)
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)
data.table::fwrite(out, file.path(DIR_OUT, "13_coloc_eqtl_igan_R.tsv"), sep = "\t")

cat("[03b] comparisons:", nrow(out), " computable:", sum(!is.na(out$PPH4)), "\n")
cat("[03b] max PPH4 =", round(max(out$PPH4, na.rm = TRUE), 4), "\n")
