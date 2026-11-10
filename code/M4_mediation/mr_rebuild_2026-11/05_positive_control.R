## =====================================================================
## 05_positive_control.R  -  pipeline positive control
## mQTL instruments -> C1GALT1 mRNA expression (eQTL Catalogue, 5 contexts)
## A valid instrument set for methylation must also carry a large, consistent
## effect on the transcript of the same gene; this calibrates sensitivity.
## out: out/10_positive_control.tsv
## =====================================================================
.this_dir <- dirname(normalizePath(sub("^--file=", "",
  commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]),
  winslash = "/", mustWork = FALSE))
source(file.path(.this_dir, "00_config.R"))

CTX <- data.frame(
  qtd   = c("QTD000356", "QTD000261", "QTD000608", "QTD000607", "QTD000606"),
  label = c("GTEx v10 whole blood", "GTEx v10 kidney cortex",
            "OneK1K B naive", "OneK1K B memory", "OneK1K B intermediate"),
  n     = c(853, 106, 876, 726, 618),
  stringsAsFactors = FALSE)

GENE <- "ENSG00000106392"   # C1GALT1

read_eqtl <- function(qtd) {
  f <- file.path(EQTL_D, sprintf("%s.7.tsv", qtd))
  d <- fread(f, sep = "\t", header = FALSE, data.table = FALSE,
             colClasses = list(character = c(1, 4, 5, 6, 19)))
  names(d)[1:19] <- c("gene", "chr", "pos", "ref", "alt", "variant", "ac", "maf", "pval",
                      "beta", "se", "type", "an", "n2", "x", "g1", "g2", "tpm", "rsid")
  d <- d[d$gene == GENE, ]
  d$beta <- as.numeric(d$beta); d$se <- as.numeric(d$se)
  d$rsid[is.na(d$rsid) | d$rsid == "NA" | d$rsid == ""] <- NA
  d[!is.na(d$se) & d$se > 0, ]
}

eq <- do.call(rbind, lapply(CTX$qtd, function(q) { x <- read_eqtl(q); x$qtd <- q; x }))
cat(sprintf("[05] cis-eQTL records for C1GALT1: %d across %d contexts\n", nrow(eq), length(CTX$qtd)))
print(table(eq$qtd))

d <- fread(file.path(DIR_OUT, "02_harmonised.tsv"), sep = "\t", data.table = FALSE)
d <- d[d$ok, ]

rows <- list()
for (i in seq_len(nrow(CTX))) {
  q <- CTX$qtd[i]
  e <- eq[eq$qtd == q & !is.na(eq$rsid), ]
  e <- e[!duplicated(e$rsid), ]
  m <- merge(d[, c("cpg", "rsid", "beta_exposure", "se_exposure", "pos")],
             e[, c("rsid", "beta", "se", "pval")], by = "rsid")
  if (nrow(m) == 0) { cat(sprintf("[05] %s: no shared variant\n", q)); next }
  for (cpg in unique(m$cpg)) {
    x <- m[m$cpg == cpg, ]
    b   <- sum(x$beta * x$beta_exposure / x$se^2) / sum(x$beta_exposure^2 / x$se^2)
    seb <- sqrt(1 / sum(x$beta_exposure^2 / x$se^2))
    rows[[length(rows) + 1]] <- data.frame(
      context = CTX$label[i], cpg = cpg, k = nrow(x),
      snps = paste(x$rsid, collapse = ";"),
      beta = b, se = seb, p = 2 * pnorm(-abs(b / seb)),
      eqtl_p_min = min(x$pval, na.rm = TRUE), stringsAsFactors = FALSE)
  }
}
pc <- do.call(rbind, rows); rownames(pc) <- NULL
pc <- pc[order(pc$cpg, pc$context), ]
cat("\n[05] positive control (methylation -> C1GALT1 expression):\n")
print(pc, digits = 4)

## consistency across contexts: same sign, all significant?
cons <- do.call(rbind, lapply(split(pc, pc$cpg), function(x)
  data.frame(cpg = x$cpg[1], n_context = nrow(x),
             n_sig = sum(x$p < 0.05), n_negative = sum(x$beta < 0),
             beta_range = paste(signif(range(x$beta), 3), collapse = " to "))))
cat("\n[05] consistency across expression contexts:\n"); print(cons)

write.table(pc,   file.path(DIR_OUT, "10_positive_control.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(cons, file.path(DIR_OUT, "10b_positive_control_consistency.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
