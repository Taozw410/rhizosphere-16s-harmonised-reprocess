#!/usr/bin/env Rscript
# 2026-10-10 Phase5b variant: NO --discard-untrimmed. Purpose: keep every sample of a MIXED-library
# project (some samples carry primers, others were already primer-free) instead of dropping the
# primer-free samples wholesale. Primers are still trimmed where present.
# 2026-10-08: cutadapt now uses --discard-untrimmed (primer-first cohort): reads without the
# primer are dropped instead of mixed in, so the ASV length frame stays uniform. truncLen values
# in the primer tables are declared in POST-trim coordinates accordingly.
# DADA2 并行管线（Windows 可用）：PSOCK 集群并行跑 filterAndTrim / dada / mergePairs
# 用法:
#   conda run -p D:\bioenv\envs\riz16s Rscript run_dada2_parallel.R <project|ALL> <raw_root> <out_root> [ref_fasta] [ncore] [nbases] [n_learn]
#
# 为什么不用 multithread=TRUE：dada2 在 Windows 上基于 fork，会直接退化为单线程
# （日志会打印 "Multithreading has been DISABLED, as forking is not supported"）。
# PSOCK 集群走 socket，Windows 可用。

suppressPackageStartupMessages({
  library(dada2)
  library(Biostrings)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("need: <project|ALL> <raw_root> <out_root> [ref_fasta] [ncore] [nbases] [n_learn]")

project   <- args[1]
raw_root  <- args[2]
out_root  <- args[3]
ref_fasta <- if (length(args) >= 4) args[4] else "D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz"
ncore     <- if (length(args) >= 5) as.integer(args[5]) else 3L
nbases    <- if (length(args) >= 6) as.numeric(args[6]) else 1e7
n_learn   <- if (length(args) >= 7) as.integer(args[7]) else 6L
primer_csv <- if (length(args) >= 8) args[8] else "D:/CodexWork/mic-re/meta/primer_table.csv"
cutadapt_exe <- "D:/bioenv/envs/riz16s/Scripts/cutadapt.exe"

primer_tab <- if (file.exists(primer_csv)) {
  # 注意：Windows 上写出的 CSV 常带 UTF-8 BOM，会把第一个列名污染成 "\ufeffproject"
  pt <- tryCatch(
    read.csv(primer_csv, stringsAsFactors = FALSE, check.names = FALSE,
             fileEncoding = "UTF-8-BOM"),
    error = function(e) read.csv(primer_csv, stringsAsFactors = FALSE, check.names = FALSE))
  names(pt) <- sub("^\ufeff", "", names(pt))
  if (!"project" %in% names(pt)) stop("primer table missing 'project' column: ", paste(names(pt), collapse = ","))
  pt
} else data.frame(project = character(0), primer_f = character(0), primer_r = character(0))

primers_for <- function(proj) {
  r <- primer_tab[primer_tab$project == proj, , drop = FALSE]
  geti <- function(col, default = 0L) {
    if (!col %in% names(r)) return(as.integer(default))
    v <- r[[col]][1]
    if (is.na(v) || !nzchar(trimws(as.character(v)))) return(as.integer(default))
    suppressWarnings(as.integer(v))
  }
  if (nrow(r) == 0) return(list(f = NA_character_, r = NA_character_, tf = 0L, tr = 0L))
  list(f = r$primer_f[1], r = r$primer_r[1],
       tf = geti("truncLen_f", 0L), tr = geti("truncLen_r", 0L))
}

# 用 cutadapt 去 5' 引物；返回 R1 带引物的百分比
run_cutadapt <- function(fF, fR, oF, oR, pF, pR) {
  # 2026-10-09: allow trimming only ONE end (R2-only projects exist: the submitter already removed the R1
  # primer, and leaving R2's primer on adds ~20-50 bp of non-biological 3' sequence to every merged ASV,
  # which depresses genus assignment). Returns the HIGHEST carriage among the ends a primer was given for,
  # so the >=50% gate works for R2-only libraries too.
  hasF <- !is.na(pF) && nzchar(pF); hasR <- !is.na(pR) && nzchar(pR)
  if ((!hasF && !hasR) || !file.exists(cutadapt_exe)) return(NA_real_)
  if (is.null(fR)) {
    if (!hasF) return(NA_real_)
    a <- c("-g", paste0("^", pF), "-o", oF, fF)
  } else {
    a <- character(0)
    if (hasF) a <- c(a, "-g", paste0("^", pF))
    if (hasR) a <- c(a, "-G", paste0("^", pR))
    a <- c(a, "-o", oF, "-p", oR, fF, fR)
  }
  out <- tryCatch(system2(cutadapt_exe, args = a, stdout = TRUE, stderr = TRUE),
                  error = function(e) NULL)
  if (is.null(out)) return(NA_real_)
  txt <- paste(out, collapse = "\n")
  pcts <- as.numeric(sub(".*\\(([0-9.]+)%\\).*", "\\1",
                         regmatches(txt, gregexpr("[0-9,]+ \\([0-9.]+%\\)", txt))[[1]]))
  if (length(pcts) == 0) return(NA_real_)
  max(pcts, na.rm = TRUE)
}

projects <- if (toupper(project) == "ALL") {
  basename(list.dirs(raw_root, recursive = FALSE))
} else {
  strsplit(project, ",")[[1]]
}

run_one <- function(proj) {
  raw_dir  <- file.path(raw_root, proj)
  out_dir  <- file.path(out_root, proj)
  filt_dir <- file.path(out_dir, "filtered")
  dir.create(filt_dir, recursive = TRUE, showWarnings = FALSE)

  logf <- file.path(out_dir, "pipeline.log")
  lc <- file(logf, open = "at")
  on.exit(try(close(lc), silent = TRUE), add = TRUE)
  say <- function(...) { m <- paste0(...); cat(m, "\n"); writeLines(m, lc); flush(lc) }
  say("===== parallel pipeline: ", proj, " =====")
  say("start: ", format(Sys.time()), "  cores=", ncore, "  nbases=", nbases)

  # 项目级跳过：已经出过完整结果的直接跳过（中断续跑时非常省时）
  if (file.exists(file.path(out_dir, "asv_taxonomy.csv")) &&
      file.exists(file.path(out_dir, "project_qc.csv"))) {
    say("SKIP (already complete)")
    return(invisible(NULL))
  }

  fns <- sort(list.files(raw_dir, pattern = "\\.fastq\\.gz$", full.names = TRUE))
  if (length(fns) == 0) { say("no fastq"); close(lc); return(invisible(NULL)) }

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
  say("mode=", ifelse(paired, "paired", "single"), "  n=", length(fnFs))

  # ---------------- 质量编码筛查 ----------------
  # 公共数据中存在质量值被替换成常数占位符的提交（所有位置、所有 reads 同一个字符），
  # 这类数据 learnErrors 必然失败（报错信息会误导成 "very few reads"）。
  quality_diversity <- function(f, nreads = 200) {
    con <- gzfile(f, "rt"); on.exit(close(con))
    ln <- readLines(con, n = nreads * 4, warn = FALSE)
    q <- ln[seq(4, length(ln), by = 4)]
    length(unique(strsplit(paste(q, collapse = ""), "")[[1]]))
  }
  qdiv <- tryCatch(quality_diversity(fnFs[1]), error = function(e) NA_integer_)
  placeholder <- !is.na(qdiv) && qdiv <= 2
  say("quality screen: ", qdiv, " distinct quality characters (first 200 reads)",
      if (placeholder) "  => PLACEHOLDER QUALITY" else "")

  fnFs_filt <- file.path(filt_dir, paste0(sample_names, "_F_filt.fastq.gz"))
  fnRs_filt <- if (paired) file.path(filt_dir, paste0(sample_names, "_R_filt.fastq.gz")) else NULL

  cl <- makeCluster(ncore, outfile = "")
  on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
  clusterEvalQ(cl, suppressPackageStartupMessages(library(dada2)))

  # ---------------- cutadapt 去引物（探针决定是否处理） ----------------
  trim_dir <- file.path(out_dir, "trimmed")
  dir.create(trim_dir, recursive = TRUE, showWarnings = FALSE)
  pr <- primers_for(proj)
  in_F <- fnFs
  in_R <- if (paired) fnRs else NULL
  trim_pct <- NA_real_
  trimmed_applied <- FALSE
  trim_csv <- file.path(out_dir, "cutadapt_report.csv")

  if (((!is.na(pr$f) && nzchar(pr$f)) || (!is.na(pr$r) && nzchar(pr$r))) && file.exists(cutadapt_exe)) {
    probe_F <- file.path(trim_dir, basename(fnFs[1]))
    probe_R <- if (paired) file.path(trim_dir, basename(fnRs[1])) else NULL
    trim_pct <- run_cutadapt(fnFs[1], if (paired) fnRs[1] else NULL,
                             probe_F, probe_R, pr$f, pr$r)
    say("cutadapt probe (sample 1): ", trim_pct, "% of R1 carry the 5' primer")
    if (!is.na(trim_pct) && trim_pct >= 50) {
      trimmed_applied <- TRUE
      say("-> primers present, trimming ALL samples")
      t0t <- Sys.time()
      tr_res <- parLapply(cl, seq_along(fnFs), function(i, fF, fR, oF, oR, pF, pR, exe) {
        hasF <- !is.na(pF) && nzchar(pF); hasR <- !is.na(pR) && nzchar(pR)
        if (is.null(fR)) a <- c("-g", paste0("^", pF), "-o", oF[i], fF[i])
        else {
          a <- character(0)
          if (hasF) a <- c(a, "-g", paste0("^", pF))
          if (hasR) a <- c(a, "-G", paste0("^", pR))
          a <- c(a, "-o", oF[i], "-p", oR[i], fF[i], fR[i])
        }
        o <- tryCatch(system2(exe, args = a, stdout = TRUE, stderr = TRUE), error = function(e) "")
        txt <- paste(o, collapse = "\n")
        mm <- regmatches(txt, regexpr("Reads? (1 )?with adapters?: *[0-9,]+ \\([0-9.]+%\\)", txt))
        data.frame(sample = basename(fF[i]),
                   pct = if (length(mm)) as.numeric(sub(".*\\(([0-9.]+)%\\).*", "\\1", mm[1])) else NA_real_)
      }, fF = fnFs, fR = if (paired) fnRs else NULL,
         oF = file.path(trim_dir, basename(fnFs)),
         oR = if (paired) file.path(trim_dir, basename(fnRs)) else NULL,
         pF = pr$f, pR = pr$r, exe = cutadapt_exe)
      tr_df <- do.call(rbind, tr_res)
      write.csv(tr_df, trim_csv, row.names = FALSE)
      say("cutadapt done in ", round(difftime(Sys.time(), t0t, units = "mins"), 1),
          " min; mean primer carriage ",
          round(mean(tr_df$pct, na.rm = TRUE), 1), "%")
      in_F <- file.path(trim_dir, basename(fnFs))
      in_R <- if (paired) file.path(trim_dir, basename(fnRs)) else NULL
    } else {
      say("-> primer not detected in sample 1; using ORIGINAL reads (no trimming)")
    }
  } else {
    say("no primer recorded / cutadapt missing; using ORIGINAL reads")
  }

  # 按项目设定的截断长度（读长超过扩增子长度时必须截断，否则配对无法合并）
  tl <- if (paired) c(pr$tf, pr$tr) else pr$tf
  say("truncLen = ", paste(tl, collapse = ","), if (all(tl == 0)) "  (no truncation)" else "  (tuned for this project)")

  # ---------------- 并行：filterAndTrim ----------------
  t0 <- Sys.time()

  filt_done <- !trimmed_applied &&
    all(file.exists(fnFs_filt)) && (!paired || all(file.exists(fnRs_filt)))
  if (filt_done && file.exists(file.path(out_dir, "filter_stats.csv"))) {
    say("filter: reuse existing filtered files")
    fout <- as.matrix(read.csv(file.path(out_dir, "filter_stats.csv"), row.names = 1))
  } else if (paired) {
    fout_list <- parLapply(cl, seq_along(fnFs), function(i, fF, fR, oF, oR, tlen) {
      dada2::filterAndTrim(fF[i], oF[i], fR[i], oR[i], truncLen = tlen, truncQ = 2,
                           maxN = 0, maxEE = c(2, 2), rm.phix = TRUE, minLen = 100,
                           compress = TRUE, multithread = FALSE, verbose = FALSE)
    }, fF = in_F, fR = in_R, oF = fnFs_filt, oR = fnRs_filt, tlen = tl)
    fout <- do.call(rbind, fout_list)
    rownames(fout) <- basename(fnFs)
    write.csv(fout, file.path(out_dir, "filter_stats.csv"))
  } else {
    fout_list <- parLapply(cl, seq_along(fnFs), function(i, fF, oF, tlen) {
      dada2::filterAndTrim(fF[i], oF[i], truncLen = tlen, truncQ = 2,
                           maxN = 0, maxEE = 2, rm.phix = TRUE, minLen = 100,
                           compress = TRUE, multithread = FALSE, verbose = FALSE)
    }, fF = in_F, oF = fnFs_filt, tlen = tl)
    fout <- do.call(rbind, fout_list)
    rownames(fout) <- basename(fnFs)
    write.csv(fout, file.path(out_dir, "filter_stats.csv"))
  }
  say("filter stage done in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")

  keep <- fout[, "reads.out"] > 0
  fnFs_filt <- fnFs_filt[keep]
  fnRs_filt <- if (paired) fnRs_filt[keep] else NULL
  sample_names <- sample_names[keep]
  say("samples kept: ", length(fnFs_filt), " / ", length(keep))

  # ---------------- 学错误率（一次性，取前 n_learn 个样本） ----------------
  errF_rds <- file.path(out_dir, "errF.rds")
  errR_rds <- file.path(out_dir, "errR.rds")
  shared_F <- file.path(out_root, "_shared", "errF.rds")
  shared_R <- file.path(out_root, "_shared", "errR.rds")
  set.seed(42)
  learn_n <- min(n_learn, length(fnFs_filt))
  t0 <- Sys.time()
  if (file.exists(errF_rds)) {
    errF <- readRDS(errF_rds); say("reuse cached errF")
  } else if (placeholder && file.exists(shared_F)) {
    errF <- readRDS(shared_F); saveRDS(errF, errF_rds)
    say("PLACEHOLDER quality -> using SHARED errF (borrowed, must be declared in Methods)")
  } else {
    say("learnErrors (F) on ", learn_n, " samples ...")
    errF <- tryCatch(
      learnErrors(fnFs_filt[seq_len(learn_n)], nbases = nbases,
                  multithread = FALSE, verbose = FALSE, MAX_CONSIST = 8),
      error = function(e) {
        say("learnErrors(F) FAILED: ", conditionMessage(e))
        if (file.exists(shared_F)) { say("-> fallback to SHARED errF"); readRDS(shared_F) } else NULL
      })
    if (is.null(errF)) { say("no usable errF for ", proj, " -- SKIP"); close(lc); return(invisible(NULL)) }
    saveRDS(errF, errF_rds)
  }
  if (paired) {
    if (file.exists(errR_rds)) {
      errR <- readRDS(errR_rds); say("reuse cached errR")
    } else if (placeholder && file.exists(shared_R)) {
      errR <- readRDS(shared_R); saveRDS(errR, errR_rds)
      say("PLACEHOLDER quality -> using SHARED errR (borrowed)")
    } else {
      say("learnErrors (R) on ", learn_n, " samples ...")
      errR <- tryCatch(
        learnErrors(fnRs_filt[seq_len(learn_n)], nbases = nbases,
                    multithread = FALSE, verbose = FALSE, MAX_CONSIST = 8),
        error = function(e) {
          say("learnErrors(R) FAILED: ", conditionMessage(e))
          if (file.exists(shared_R)) { say("-> fallback to SHARED errR"); readRDS(shared_R) } else NULL
        })
      if (is.null(errR)) { say("no usable errR for ", proj, " -- SKIP"); close(lc); return(invisible(NULL)) }
      saveRDS(errR, errR_rds)
    }
  } else errR <- NULL
  say("error learning done in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")

  # ---------------- 并行：dada + mergePairs ----------------
  t0 <- Sys.time()
  deno_rds <- file.path(out_dir, "deno_list.rds")
  if (file.exists(deno_rds) && !trimmed_applied) {
    say("denoise: reuse cached deno_list.rds")
    deno_list <- readRDS(deno_rds)
  } else if (paired) {
    deno_list <- parLapply(cl, seq_along(fnFs_filt), function(i, fF, fR, eF, eR) {
      dF <- dada2::dada(fF[i], err = eF, multithread = FALSE, pool = FALSE, verbose = FALSE)
      dR <- dada2::dada(fR[i], err = eR, multithread = FALSE, pool = FALSE, verbose = FALSE)
      dada2::mergePairs(dF, fF[i], dR, fR[i], verbose = FALSE)
    }, fF = fnFs_filt, fR = fnRs_filt, eF = errF, eR = errR)
    names(deno_list) <- sample_names
    saveRDS(deno_list, deno_rds)
  } else {
    deno_list <- parLapply(cl, seq_along(fnFs_filt), function(i, fF, eF) {
      dada2::dada(fF[i], err = eF, multithread = FALSE, pool = FALSE, verbose = FALSE)
    }, fF = fnFs_filt, eF = errF)
    names(deno_list) <- sample_names
    saveRDS(deno_list, deno_rds)
  }
  say("denoise+merge done in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")

  # ---------------- 串行：建表 / 去嵌合 / 分类学 ----------------
  seqtab <- makeSequenceTable(deno_list)
  say("seqtab: ", paste(dim(seqtab), collapse = " x "))
  seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                      multithread = FALSE, verbose = FALSE)
  say("after chimera: ", paste(dim(seqtab_nochim), collapse = " x "),
      "  reads kept=", round(sum(seqtab_nochim) / sum(seqtab), 4))

  write.csv(t(seqtab_nochim), file.path(out_dir, "asv_table.csv"))
  saveRDS(seqtab_nochim, file.path(out_dir, "seqtab_nochim.rds"))

  if (paired) {
    track <- data.frame(sample = sample_names,
                        input = fout[keep, "reads.in"],
                        filtered = fout[keep, "reads.out"],
                        merged = sapply(deno_list, function(x) sum(getUniques(x))),
                        nonchim = rowSums(seqtab_nochim), row.names = NULL)
  } else {
    track <- data.frame(sample = sample_names,
                        input = fout[keep, "reads.in"],
                        filtered = fout[keep, "reads.out"],
                        denoised = sapply(deno_list, function(x) sum(getUniques(x))),
                        nonchim = rowSums(seqtab_nochim), row.names = NULL)
  }
  write.csv(track, file.path(out_dir, "read_tracking.csv"), row.names = FALSE)

  taxa <- tryCatch(
    assignTaxonomy(seqtab_nochim, ref_fasta, multithread = FALSE, verbose = FALSE),
    error = function(e) { say("assignTaxonomy FAILED: ", conditionMessage(e)); NULL })
  if (!is.null(taxa)) {
    write.csv(taxa, file.path(out_dir, "asv_taxonomy.csv"))
    say("taxonomy: genus assigned ", sum(!is.na(taxa[, "Genus"])), " / ", nrow(taxa))
  }
  write.csv(data.frame(project = proj, mode = ifelse(paired, "paired", "single"),
                       n_samples = length(fnFs_filt), quality_chars = qdiv,
                       placeholder_quality = placeholder,
                       shared_err_matrix = placeholder,
                       n_asv = ncol(seqtab_nochim), row.names = NULL),
            file.path(out_dir, "project_qc.csv"), row.names = FALSE)
  say("DONE ", proj, " at ", format(Sys.time()))
  close(lc)
  invisible(NULL)
}

cat("projects:", paste(projects, collapse = ", "), "\n")
for (p in projects) {
  cat("\n########## ", p, " ##########\n")
  try(run_one(p))
}
cat("\nall done:", format(Sys.time()), "\n")
