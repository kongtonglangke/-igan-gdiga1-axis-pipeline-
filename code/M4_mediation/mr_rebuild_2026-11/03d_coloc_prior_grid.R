## =====================================================================
## 03d_coloc_prior_grid.R  -  prior-sensitivity grid for the mQTL x eQTL call
##
## Re-runs the cg19473623 x kidney-cortex comparison over
##   shared prior effect variance W = 0.01, 0.02, 0.04, 0.10
##   p12 = 1e-6, 1e-5, 1e-4
## with coloc::coloc.abf, so the stability statement in the Results is backed by
## the same implementation as the main analysis.
##
## W is set through the expression-side sd.prior * sdY product (coloc fixes
## sd.prior at 0.15 for a quantitative trait): sdY = sqrt(W) / 0.15.
##
## out: out/15_coloc_prior_grid.tsv
## =====================================================================
.this_dir <- dirname(normalizePath(sub("^--file=", "",
  commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]),
  winslash = "/", mustWork = FALSE))
source(file.path(.this_dir, "00_config.R"))

suppressPackageStartupMessages(library(coloc))

CPG <- "cg19473623"
QTD <- "QTD000261"          # GTEx v10 kidney cortex
ENSG <- "ENSG00000106392"

mq <- data.table::fread(file.path(MQLT_D, sprintf("cpg_full_assoc_%s.tsv", CPG)),
                        sep = "\t", header = TRUE, showProgress = FALSE)
mq <- mq[, .(rsid, beta = as.numeric(beta_are_a1), se = as.numeric(se),
             a1 = toupper(a1), a2 = toupper(a2))]
mq <- mq[!is.na(rsid) & rsid != "" & !is.na(beta) & !is.na(se) & se > 0]
mq <- mq[!duplicated(rsid)]

eq <- data.table::fread(file.path(ROOT, "data/eqtl_for_coloc", sprintf("%s.7.tsv", QTD)),
                        sep = "\t", header = TRUE, showProgress = FALSE)
eq <- eq[gene == ENSG]
eq[, `:=`(alt = toupper(alt), ref = toupper(ref))]
eq <- eq[!duplicated(rsid)]

m <- merge(mq, eq, by = "rsid", suffixes = c(".m", ".e"))
same <- m$a1 == m$alt & m$a2 == m$ref
flip <- m$a1 == m$ref & m$a2 == m$alt
keep <- same | flip
m <- m[keep, , drop = FALSE]
b <- as.numeric(m[["beta.e"]])
f <- as.logical(flip[keep])
b[f] <- -b[f]
m[["beta.e"]] <- b

rows <- list()
for (W in c(0.01, 0.02, 0.04, 0.10)) {
  sdY <- sqrt(W) / 0.15
  for (p12 in c(1e-6, 1e-5, 1e-4)) {
    d1 <- list(snp = m$rsid, beta = m[["beta.e"]], varbeta = as.numeric(m[["se.e"]])^2,
               position = m$pos, type = "quant", N = 700, sdY = sdY)
    d2 <- list(snp = m$rsid, beta = as.numeric(m[["beta.m"]]),
               varbeta = as.numeric(m[["se.m"]])^2,
               position = m$pos, type = "quant", N = 25000, sdY = 1)
    r <- coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = p12)
    s <- r$summary
    rows[[length(rows) + 1]] <- data.frame(
      W = W, p12 = p12, n_snp = as.integer(s[["nsnps"]]),
      PPH3 = s[["PP.H3.abf"]], PPH4 = s[["PP.H4.abf"]])
  }
}

out <- do.call(rbind, rows)
data.table::fwrite(out, file.path(DIR_OUT, "15_coloc_prior_grid.tsv"), sep = "\t")
cat("[03d] cg19473623 x kidney cortex,", nrow(out), "prior combinations\n")
cat("[03d] PPH4 range:", sprintf("%.3f - %.4f", min(out$PPH4), max(out$PPH4)), "\n")
print(out)
