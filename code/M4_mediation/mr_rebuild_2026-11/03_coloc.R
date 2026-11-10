## =====================================================================
## 03_coloc.R  -  colocalisation mQTL x IgAN at the C1GALT1 locus
## priors: p1 = p2 = 1e-4, p12 = 1e-5 (as pre-specified in the study protocol)
## out: out/06_coloc_mqtl_igan.tsv   out/06b_coloc_sensitivity.tsv
## =====================================================================
.this_dir <- dirname(normalizePath(sub("^--file=", "",
  commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]),
  winslash = "/", mustWork = FALSE))
source(file.path(.this_dir, "00_config.R"))

mqtl <- fread(file.path(DIR_OUT, "01_instrument_candidates.tsv"), sep = "\t", data.table = FALSE)
igan <- fread(file.path(DIR_DATA, "igan_chr7_region.tsv"), sep = "\t", header = TRUE, data.table = FALSE)
names(igan) <- c("chr", "pos", "ea", "oa", "beta", "se", "eaf", "pval", "vid", "hm_conv", "hm_code", "rsid")
igan <- igan[!is.na(igan$rsid) & !is.na(igan$beta) & !is.na(igan$se) & igan$se > 0, ]
igan <- igan[!duplicated(igan$rsid), ]

S_CC <- N_IGAN_CASE / N_IGAN          # case proportion for the cc likelihood

res <- list(); sens <- list()
for (cpg in CPG_INFO$cpg) {
  x <- mqtl[mqtl$cpg == cpg & !is.na(mqtl$rsid), ]
  x <- x[!duplicated(x$rsid), ]
  m <- merge(x[, c("rsid", "pos", "a1", "a2", "beta", "se", "maf", "eaf", "samplesize")],
             igan[, c("rsid", "ea", "oa", "beta", "se", "eaf")], by = "rsid")
  if (nrow(m) < 10) { cat(sprintf("[03] %s: only %d shared SNPs - skipped\n", cpg, nrow(m))); next }

  ## allele alignment to the mQTL effect allele
  ok <- toupper(m$a1) == toupper(m$ea) & toupper(m$a2) == toupper(m$oa)
  flip <- toupper(m$a1) == toupper(m$oa) & toupper(m$a2) == toupper(m$ea)
  m <- m[ok | flip, ]
  m$by <- ifelse(toupper(m$a1) == toupper(m$ea), m$beta.y, -m$beta.y)
  m$byse <- m$se.y
  m$bx <- m$beta.x; m$bxse <- m$se.x
  ## palindromic SNPs are uninformative for allele alignment; drop them when MAF is high
  palin <- (toupper(m$a1) == "A" & toupper(m$a2) == "T") | (toupper(m$a1) == "T" & toupper(m$a2) == "A") |
           (toupper(m$a1) == "C" & toupper(m$a2) == "G") | (toupper(m$a1) == "G" & toupper(m$a2) == "C")
  m <- m[!(palin & pmin(m$eaf.x, 1 - m$eaf.x) > 0.42), ]

  d1 <- list(beta = m$bx, varbeta = m$bxse^2, snp = m$rsid, position = m$pos %% 1e6,
             type = "quant", N = median(m$samplesize, na.rm = TRUE), MAF = m$maf)
  d2 <- list(beta = m$by, varbeta = m$byse^2, snp = m$rsid, position = m$pos %% 1e6,
             type = "cc", N = N_IGAN, s = S_CC, MAF = m$eaf.y)
  cc <- suppressMessages(coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5))
  s  <- cc$summary
  res[[cpg]] <- data.frame(cpg = cpg, n_snp = length(m$rsid),
    PPH0 = s[["PP.H0.abf"]], PPH1 = s[["PP.H1.abf"]], PPH2 = s[["PP.H2.abf"]],
    PPH3 = s[["PP.H3.abf"]], PPH4 = s[["PP.H4.abf"]],
    top_snp = cc$results$snp[which.max(cc$results[[grep("H4", names(cc$results), value = TRUE)[1]]])],
    stringsAsFactors = FALSE)

  ## sensitivity: alternative prior p12 = 1e-4 and 1e-6
  for (p12 in c(1e-6, 1e-4)) {
    cc2 <- suppressMessages(coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = p12))
    sens[[paste(cpg, p12)]] <- data.frame(cpg = cpg, p12 = p12,
      PPH3 = cc2$summary[["PP.H3.abf"]], PPH4 = cc2$summary[["PP.H4.abf"]])
  }
  cat(sprintf("[03] %s: nSNP=%d  PPH1=%.4f  PPH3=%.4f  PPH4=%.4f\n", cpg, length(m$rsid),
              s[["PP.H1.abf"]], s[["PP.H3.abf"]], s[["PP.H4.abf"]]))
}
coloc_res <- do.call(rbind, res)
sens_res  <- do.call(rbind, sens)

## CpGs with no cis-mQTL record
miss <- setdiff(CPG_INFO$cpg, coloc_res$cpg)
if (length(miss)) coloc_res <- rbind(coloc_res, data.frame(cpg = miss, n_snp = 0,
  PPH0 = NA, PPH1 = NA, PPH2 = NA, PPH3 = NA, PPH4 = NA, top_snp = NA))

write.table(coloc_res, file.path(DIR_OUT, "06_coloc_mqtl_igan.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(sens_res,  file.path(DIR_OUT, "06b_coloc_sensitivity.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("\n[03] coloc summary:\n"); print(coloc_res, digits = 3)
cat(sprintf("\n[03] max PPH4 = %.4f ; min PPH3 among computable = %.4f\n",
            max(coloc_res$PPH4, na.rm = TRUE), min(coloc_res$PPH3, na.rm = TRUE)))
