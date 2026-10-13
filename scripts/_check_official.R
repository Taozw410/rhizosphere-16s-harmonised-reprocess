suppressPackageStartupMessages(library(DECIPHER)); suppressPackageStartupMessages(library(Biostrings))
p <- "D:/CodexWork/mic-re/ref_region/SILVA_SSU_r138_2_2024.RData"
t0 <- Sys.time()
e <- new.env()
objs <- load(p, envir = e)
cat("loaded objects:", paste(objs, collapse = ", "), "\n")
cat("load time (min):", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
for (o in objs) {
  x <- get(o, envir = e)
  cat(sprintf("  %-16s class=%s  length=%s  size=%.1f MB\n", o, paste(class(x), collapse = "/"),
              ifelse(is.null(length(x)), "-", length(x)), as.numeric(object.size(x)) / 1024^2))
}
ts <- get(objs[1], envir = e)
if (inherits(ts, "Train")) cat("taxa groups:", length(ts$taxa), "\n")
cat("mem used (GB):", round(sum(gc()[, 2]), 2), "\n")
# 用 3 条真实 ASV 试 IdTaxa
rd <- read.csv("D:/CodexWork/mic-re/out_stage2/PRJNA633917/asv_table.csv", row.names = 1, check.names = FALSE)
sq <- DNAStringSet(rownames(rd)[1:3]); names(sq) <- paste0("q", 1:3)
ids <- IdTaxa(sq, ts, strand = "top", threshold = 60, processors = 1, verbose = FALSE)
for (i in seq_along(ids)) cat(sprintf("  %s -> %s\n", names(sq)[i], paste(ids[[i]]$taxon, collapse = "; ")))
