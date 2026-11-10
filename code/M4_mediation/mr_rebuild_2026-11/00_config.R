## =====================================================================
## 00_config.R  -  global settings (sourced by every script)
## Project: C1GALT1 cis-mQTL -> IgAN two-sample Mendelian randomisation
## Rebuilt 2026-11-10.  ASCII-only paths on purpose (Windows + Chinese
## path breaks Rscript source() decoding).
##
## Paths are resolved at run time, never hard-coded:
##   MR_ROOT     analysis root; inputs go to <MR_ROOT>/data, tables to <MR_ROOT>/out
##               (default: the directory that contains this R/ tree)
##   MQTL_DIR    the five GoDMC cis-mQTL tables `cpg_full_assoc_<CpG>.tsv`
##               (default: <MR_ROOT>/data/mqtl)
##   EQTL_CACHE  directory holding the eQTL Catalogue region cache
##               (default: <MR_ROOT>/data/eQTL_region_cache)
##   MR_FIG_DIR  figure output directory (default: <MR_ROOT>/fig)
##   MR_DATA_DIR directory holding the two IgAN input extracts
##               (default: <MR_ROOT>/data)
## =====================================================================

.mr_script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) == 1L) dirname(normalizePath(f, winslash = "/", mustWork = FALSE)) else "R"
}

ROOT   <- Sys.getenv("MR_ROOT", unset = file.path(dirname(.mr_script_dir()), "."))
ROOT   <- normalizePath(ROOT, winslash = "/", mustWork = FALSE)
MQLT_D <- Sys.getenv("MQTL_DIR", unset = file.path(ROOT, "data/mqtl"))
EQTL_D <- Sys.getenv("EQTL_CACHE", unset = file.path(ROOT, "data/eQTL_region_cache"))
FIG_D  <- Sys.getenv("MR_FIG_DIR", unset = file.path(ROOT, "fig"))

DIR_DATA <- Sys.getenv("MR_DATA_DIR", unset = file.path(ROOT, "data"))
DIR_OUT <- file.path(ROOT, "out"); DIR_FIG <- FIG_D; DIR_LOG <- file.path(ROOT, "log")
for (d in c(DIR_DATA, DIR_OUT, DIR_FIG, DIR_LOG))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

## ---- analysis parameters (explicit, for reproducibility) ----
P_INSTR  <- 5e-8
P_RELAX  <- 1e-5
F_MIN    <- 10
CLUMP_R2 <- 0.001
CLUMP_KB <- 10000

## ---- sample sizes (from source metadata; do not edit) ----
N_IGAN_CASE <- 15658
N_IGAN_CTRL <- 653143 - 15658
N_IGAN      <- 653143
N_IGAN_EFF  <- 4 / (1/N_IGAN_CASE + 1/N_IGAN_CTRL)

## ---- CpG annotation ----
CPG_INFO <- data.frame(
  cpg   = c("cg19603390", "cg19473623", "cg17994788", "cg04827551", "cg16101574"),
  pos37 = c(7222222, 7224869, 7261594, 7268805, 7291514),
  gene  = "C1GALT1",
  relation = c("intragenic", "intragenic", "intragenic", "intragenic", "3'downstream<5kb"),
  stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2); library(data.table); library(dplyr); library(tidyr)
  library(patchwork); library(MendelianRandomization); library(MRPRESSO)
  library(RadialMR); library(coloc); library(scales)
})

theme_sci <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(axis.line  = element_line(linewidth = 0.4, colour = "black"),
          axis.ticks = element_line(linewidth = 0.4, colour = "black"),
          axis.text  = element_text(colour = "black", size = base_size),
          axis.title = element_text(colour = "black", size = base_size + 0.5),
          plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0),
          plot.subtitle = element_text(size = base_size - 0.5, colour = "grey25"),
          strip.background = element_rect(fill = "grey95", colour = NA),
          strip.text = element_text(face = "bold", size = base_size),
          legend.key.size = unit(3.2, "mm"),
          legend.text = element_text(size = base_size - 0.5),
          legend.title = element_text(size = base_size - 0.5, face = "bold"),
          plot.tag = element_text(face = "bold", size = 9),
          plot.margin = margin(4, 6, 4, 4))
}
theme_set(theme_sci())

PAL <- c(mQTL = "#B2182B", IgAN = "#2166AC", eQTL = "#1B7837",
         neut = "#4D4D4D", warn = "#E08214", grey = "#BDBDBD")

pct  <- function(x, d = 1) paste0(formatC(100 * x, format = "f", digits = d), "%")
fmtp <- function(p) ifelse(p < 1e-4, formatC(p, format = "e", digits = 2),
                           formatC(p, format = "f", digits = 4))
