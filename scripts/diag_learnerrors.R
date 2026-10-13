suppressPackageStartupMessages(library(dada2))
f <- "D:/CodexWork/mic-re/out/filtered/PRJNA1049404/SRR27115409_1_filt.fastq.gz"
r <- "D:/CodexWork/mic-re/out/filtered/PRJNA1049404/SRR27115409_2_filt.fastq.gz"
cat("file exists:", file.exists(f), file.exists(r), "\n")
err <- tryCatch(learnErrors(c(f), nbases=1e7, multithread=FALSE, verbose=1), error=function(e) {cat("ERR:", conditionMessage(e), "\n"); NULL})
if (!is.null(err)) { cat("learnErrors OK\n"); print(head(err$err_out)) }
