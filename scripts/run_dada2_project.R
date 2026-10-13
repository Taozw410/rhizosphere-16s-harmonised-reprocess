#!/usr/bin/env Rscript
# DADA2 单项目管线：过滤 -> 去噪 -> 合并 -> 去嵌合体 -> 分类学注释
# 用法:
#   conda run -p D:\bioenv\envs\riz16s Rscript run_dada2_project.R <project> <raw_root> <out_root> [ref_fasta] [nthreads] [nbases] [n_learn]
#
# 说明（2026-10-05 实测）：
#   * Windows 上 dada2 的 multithread 无效（forking 不支持），所有步骤实际单线程；
#     真正的并行只能在脚本外用多进程做。
#   * learnErrors 默认 nbases=1e8 在本机会跑到 40 分钟以上；默认改为 1e7，并只取 6 个样本学习。
#   * 错误率矩阵会缓存为 rds，重跑不再重复学习。

suppressPackageStartupMessages({
  library(dada2)
  library(Biostrings)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("need: <project> <raw_root> <out_root> [ref_fasta] [nthreads]")

project   <- args[1]
raw_root  <- args[2]
out_root  <- args[3]
ref_fasta <- if (length(args) >= 4) args[4] else "D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz"
nthreads  <- if (length(args) >= 5) as.integer(args[5]) else 2L
nbases    <- if (length(args) >= 6) as.numeric(args[6]) else 1e7
n_learn   <- if (length(args) >= 7) as.integer(args[7]) else 6L

raw_dir  <- file.path(raw_root, project)
filt_dir <- file.path(out_root, "filtered", project)
out_dir  <- file.path(out_root, project)
dir.create(filt_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)

log_con <- file(file.path(out_dir, "pipeline.log"), open = "wt")
sink(log_con, split = TRUE)
on.exit({ sink(); close(log_con) }, add = TRUE)

cat("=== DADA2 pipeline ===\n")
cat("project  :", project, "\n")
cat("raw_dir  :", raw_dir, "\n")
cat("out_dir  :", out_dir, "\n")
cat("threads  :", nthreads, "\n")
cat("nbases    :", nbases, "\n")
cat("n_learn   :", n_learn, "\n")
cat("start    :", format(Sys.time()), "\n\n")

# ---------------------------------------------------------------- 输入文件
fns <- sort(list.files(raw_dir, pattern = "\\.fastq\\.gz$", full.names = TRUE))
if (length(fns) == 0) stop("no fastq.gz found in ", raw_dir)

paired <- any(grepl("_1\\.fastq\\.gz$", fns)) && any(grepl("_2\\.fastq\\.gz$", fns))
if (paired) {
  fnFs <- sort(fns[grepl("_1\\.fastq\\.gz$", fns)])
  fnRs <- sort(fns[grepl("_2\\.fastq\\.gz$", fns)])
  sample_names <- sub("_1\\.fastq\\.gz$", "", basename(fnFs))
} else {
  fnFs <- fns
  fnRs <- NULL
  sample_names <- sub("\\.fastq\\.gz$", "", basename(fnFs))
}
cat("mode           :", ifelse(paired, "paired-end", "single-end"), "\n")
cat("n samples      :", length(fnFs), "\n\n")

# ---------------------------------------------------------------- 过滤
fnFs_filt <- file.path(filt_dir, paste0(sample_names, "_F_filt.fastq.gz"))
fnRs_filt <- if (paired) file.path(filt_dir, paste0(sample_names, "_R_filt.fastq.gz")) else NULL

cat("--- filterAndTrim ---\n")
fout <- filterAndTrim(
  fwd = fnFs, filt = fnFs_filt,
  rev = fnRs, filt.rev = fnRs_filt,
  truncLen = 0, truncQ = 2, maxN = 0, maxEE = c(2, 2),
  rm.phix = TRUE, minLen = 100, compress = TRUE,
  multithread = nthreads, verbose = TRUE
)
print(fout)
cat("\n")

keep <- fout[, "reads.out"] > 0
if (sum(keep) == 0) stop("all samples empty after filtering")
fnFs_filt <- fnFs_filt[keep]
fnRs_filt <- if (paired) fnRs_filt[keep] else NULL
sample_names <- sample_names[keep]
cat("samples kept :", length(fnFs_filt), "of", length(keep), "\n\n")

# ---------------------------------------------------------------- 去噪
cat("--- learnErrors / dada ---\n")
# 只取前 n_learn 个样本学习错误率（Windows 单线程，全量样本会极慢）
set.seed(42)
learn_n <- min(n_learn, length(fnFs_filt))
cat("learning error rates on", learn_n, "of", length(fnFs_filt), "samples\n")

errF_rds <- file.path(out_dir, "errF.rds")
errR_rds <- file.path(out_dir, "errR.rds")

if (file.exists(errF_rds)) {
  cat("reusing cached errF\n"); errF <- readRDS(errF_rds)
} else {
  errF <- learnErrors(fnFs_filt[seq_len(learn_n)], nbases = nbases,
                      multithread = FALSE, verbose = 1, MAX_CONSIST = 8)
  saveRDS(errF, errF_rds)
}
if (paired) {
  if (file.exists(errR_rds)) {
    cat("reusing cached errR\n"); errR <- readRDS(errR_rds)
  } else {
    errR <- learnErrors(fnRs_filt[seq_len(learn_n)], nbases = nbases,
                        multithread = FALSE, verbose = 1, MAX_CONSIST = 8)
    saveRDS(errR, errR_rds)
  }
}

denoF <- dada(fnFs_filt, err = errF, multithread = FALSE, pool = "pseudo", verbose = 1)
if (paired) denoR <- dada(fnRs_filt, err = errR, multithread = FALSE, pool = "pseudo", verbose = 1)

if (paired) {
  mergers <- mergePairs(denoF, fnFs_filt, denoR, fnRs_filt, verbose = TRUE)
  seqtab <- makeSequenceTable(mergers)
} else {
  seqtab <- makeSequenceTable(denoF)
}
cat("\nseqtab dim (samples x ASVs):", paste(dim(seqtab), collapse = " x "), "\n\n")

# ---------------------------------------------------------------- 去嵌合体
cat("--- removeBimeraDenovo ---\n")
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = nthreads, verbose = TRUE)
cat("after chimera removal:", paste(dim(seqtab_nochim), collapse = " x "), "\n")
cat("fraction of reads kept after chimera removal:",
    round(sum(seqtab_nochim) / sum(seqtab), 4), "\n\n")

# ---------------------------------------------------------------- 分类学
# 先把 ASV 表落盘：分类学注释依赖外部库，失败时不能让前面的工作白费
write.csv(t(seqtab_nochim), file.path(out_dir, "asv_table.csv"))
saveRDS(seqtab_nochim, file.path(out_dir, "seqtab_nochim.rds"))
cat("saved asv_table.csv and seqtab_nochim.rds\n")

cat("--- assignTaxonomy (SILVA) ---\n")
taxa <- tryCatch(
  assignTaxonomy(seqtab_nochim, ref_fasta, multithread = FALSE, verbose = TRUE),
  error = function(e) {
    cat("assignTaxonomy FAILED:", conditionMessage(e), "\n")
    NULL
  }
)
if (is.null(taxa)) {
  cat("taxonomy skipped；ASV 表已保存，可稍后单独补跑分类学\n")
  quit(status = 0)
}
cat("\n")

# ---------------------------------------------------------------- 输出
write.csv(taxa, file.path(out_dir, "asv_taxonomy.csv"))

track <- data.frame(
  sample   = sample_names,
  input    = fout[keep, "reads.in"],
  filtered = fout[keep, "reads.out"],
  row.names = NULL
)
if (paired) {
  track$denoised_F <- sapply(denoF, function(x) sum(getUniques(x)))
  track$denoised_R <- sapply(denoR, function(x) sum(getUniques(x)))
  track$merged     <- sapply(mergers, function(x) sum(getUniques(x)))
} else {
  track$denoised   <- sapply(denoF, function(x) sum(getUniques(x)))
}
track$nonchim <- rowSums(seqtab_nochim)
write.csv(track, file.path(out_dir, "read_tracking.csv"), row.names = FALSE)

cat("\n=== read tracking (head) ===\n")
print(head(track, 10))

cat("\n=== taxonomy level coverage ===\n")
for (lvl in 1:min(7, ncol(taxa))) {
  nm <- colnames(taxa)[lvl]
  cat(sprintf("  %-10s assigned: %d / %d (%.1f%%)\n", nm,
              sum(!is.na(taxa[, lvl])), nrow(taxa),
              100 * sum(!is.na(taxa[, lvl])) / nrow(taxa)))
}

cat("\n=== top 15 genera ===\n")
gen <- if ("Genus" %in% colnames(taxa)) taxa[, "Genus"] else NA
if (!all(is.na(gen))) {
  ab <- tapply(colSums(seqtab_nochim), gen, sum)
  ab <- sort(ab, decreasing = TRUE)
  print(head(round(100 * ab / sum(ab), 2), 15))
}

cat("\ndone:", format(Sys.time()), "\n")
