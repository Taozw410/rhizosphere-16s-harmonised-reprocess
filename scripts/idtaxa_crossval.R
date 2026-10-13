#!/usr/bin/env Rscript
# IDTAXA cross-validation (pre-registered section 3).  Source is ASCII-only ON PURPOSE:
# R on Windows reads this file as GBK, so Chinese text here gets mangled and breaks parsing (AGENTS 10.5).
# usage: Rscript idtaxa_crossval.R <train.fa> <out_dir> <asv_per_project> <project_list.tsv>
args <- commandArgs(trailingOnly = TRUE)
train_fa <- args[1]
out_dir  <- args[2]
n_asv    <- if (length(args) >= 3) as.integer(args[3]) else 300L
list_tsv <- if (length(args) >= 4) args[4] else "D:/CodexWork/mic-re/out50k/idtaxa/project_list.tsv"
root <- "D:/CodexWork/mic-re"
suppressPackageStartupMessages(library(DECIPHER))
suppressPackageStartupMessages(library(Biostrings))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "  ", sprintf(...), "\n", sep = "")

official <- file.path(root, "ref_region", "SILVA_SSU_r138_2_2024.RData")
tr_rds <- file.path(out_dir, "idtaxa_training.rds")
if (file.exists(official)) {
  say("loading DECIPHER official training set: %s", basename(official))
  e <- new.env(); objs <- load(official, envir = e); tr <- get(objs[1], envir = e); rm(e); gc(verbose = FALSE)
  say("official training ready: %d taxa groups", length(tr$taxa))
} else if (file.exists(tr_rds)) {
  tr <- readRDS(tr_rds); say("reuse cached self-built training object")
} else {
  seqs <- readDNAStringSet(train_fa)
  tax <- paste0("Root;", sub("^\\S+\\s+", "", names(seqs)))
  say("training sequences: %d", length(seqs))
  tr <- LearnTaxa(seqs, tax, verbose = FALSE)
  saveRDS(tr, tr_rds)
}
# tr$taxonomy holds full PATH strings ("Root;Bacteria;...;Genus;") and tr$ranks the rank of each path's
# last element, so the rank of an assignment is found by matching the joined path (with trailing ";").
genus_of <- function(tv) {
  if (length(tv) < 2) return("")
  g <- tv[length(tv)]
  if (grepl("^g__", g)) return(sub("^g__", "", g))
  i <- match(paste0(paste(tv, collapse = ";"), ";"), tr$taxonomy)
  if (is.na(i)) return("")
  r <- tr$ranks[i]
  if (!is.na(r) && r == "genus") return(g)
  ""
}
plist <- read.delim(list_tsv, stringsAsFactors = FALSE, quote = "")
say("projects to annotate: %d", nrow(plist))
res <- NULL
asv_calls <- NULL
for (i in seq_len(nrow(plist))) {
  p <- plist$project[i]; bs <- plist$base[i]
  pd <- file.path(root, bs, p)
  ans <- try({
    ab <- read.csv(file.path(pd, "asv_table.csv"), row.names = 1, check.names = FALSE)
    tx <- read.csv(file.path(pd, "asv_taxonomy.csv"), row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
    tot <- rowSums(ab, na.rm = TRUE)
    keep <- names(sort(tot, decreasing = TRUE))[1:min(n_asv, length(tot))]
    sq <- DNAStringSet(keep); names(sq) <- sprintf("a%05d", seq_along(sq))
    ids <- IdTaxa(sq, tr, strand = "top", threshold = 60, processors = NULL, verbose = FALSE)
    gid <- vapply(ids, function(x) genus_of(x$taxon), character(1))
    gdd <- vapply(keep, function(s) { v <- tx[s, "Genus"]; if (is.null(v) || is.na(v) || toupper(trimws(v)) == "NA") "" else trimws(v) }, character(1))
    both <- gid != "" & gdd != ""
    agree <- if (sum(both) > 0) 100 * sum(gid[both] == gdd[both]) / sum(both) else NA_real_
    c(n = length(keep), both = sum(both), agree = agree,
      id_rate = 100 * sum(gid != "") / length(keep), dd_rate = 100 * sum(gdd != "") / length(keep),
      read_share = 100 * sum(tot[keep]) / sum(tot))
  }, silent = TRUE)
  if (inherits(ans, "try-error")) { say("%s FAILED: %s", p, substr(as.character(ans), 1, 110)); next }
  asv_calls <- rbind(asv_calls, data.frame(project = p, asv = keep, reads = as.integer(tot[keep]),
                                           idtaxa_genus = gid, dada2_genus = gdd,
                                           agree = ifelse(gid != "" & gdd != "", as.character(gid == gdd), "")))
  res <- rbind(res, data.frame(project = p, n_asv = ans[["n"]], both = ans[["both"]],
                               agree_pct = round(ans[["agree"]], 1), idtaxa_genus_rate = round(ans[["id_rate"]], 1),
                               dada2_genus_rate = round(ans[["dd_rate"]], 1), read_share = round(ans[["read_share"]], 1)))
  say("%s  n=%d  both=%d  agree=%.1f%%  idtaxa=%.1f%%  dada2=%.1f%%", p, ans[["n"]], ans[["both"]],
      ans[["agree"]], ans[["id_rate"]], ans[["dd_rate"]])
}
if (!is.null(res)) {
  write.csv(res, file.path(out_dir, "idtaxa_vs_dada2.csv"), row.names = FALSE)
  # 逐条 ASV 的两法属级比对：核心-边缘框架的敏感性分析要据此生成"两法一致的属"清单（§3）
  if (!is.null(asv_calls)) {
    write.csv(asv_calls, file.path(out_dir, "idtaxa_asv_calls.csv"), row.names = FALSE)
    say("wrote per-ASV calls: %d rows", nrow(asv_calls))
  }
  w <- res$both
  overall <- if (sum(w) > 0) sum(res$agree_pct * w / 100) / sum(w) * 100 else NA_real_
  d_rate <- mean(res$idtaxa_genus_rate) - mean(res$dada2_genus_rate)
  say("OVERALL read-weighted genus agreement = %.1f%%  (projects=%d)", overall, nrow(res))
  say("mean genus rate: idtaxa=%.1f%%  dada2=%.1f%%  diff=%+.1fpp", mean(res$idtaxa_genus_rate), mean(res$dada2_genus_rate), d_rate)
  verdict <- if (!is.na(overall) && overall >= 80 && abs(d_rate) <= 10) "CONSISTENT (agreement >=80% and rate diff <=10pp)" else "NOT CONSISTENT -> use only genera agreed by both methods for sensitivity (section 3)"
  cat("VERDICT:", verdict, "\n")
  cat(sprintf("MEAN_READ_SHARE_OF_TOP300=%.1f%%\n", mean(res$read_share)))
}
