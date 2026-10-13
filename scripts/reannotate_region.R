#!/usr/bin/env Rscript
# 用区域匹配参考重新注释某个项目的 ASV（只重跑注释，不重跑 ASV 推断）
suppressPackageStartupMessages({ library(dada2); library(Biostrings) })
args <- commandArgs(trailingOnly = TRUE)
at <- args[1]; ref <- args[2]; out <- args[3]
con <- file(at, open = "rt")
readLines(con, n = 1)               # header
seqs <- character(0)
repeat {
  ln <- readLines(con, n = 1)
  if (length(ln) == 0) break
  if (!nzchar(ln)) next
  s <- strsplit(ln, ",", fixed = TRUE)[[1]][1]
  seqs <- c(seqs, gsub('"', '', s))
}
close(con)
proj <- basename(dirname(at))
cat(sprintf("== %s : %d ASVs ; ref=%s\n", proj, length(seqs), basename(ref)))
tax <- tryCatch(assignTaxonomy(DNAStringSet(seqs), ref, multithread = FALSE, verbose = FALSE),
                error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL })
if (!is.null(tax)) {
  write.csv(tax, out)
  g <- sum(!is.na(tax[, "Genus"]) & tax[, "Genus"] != "NA")
  cat(sprintf("%s: genus=%d / %d = %.1f%%\n", proj, g, length(seqs), 100 * g / length(seqs)))
}
