p <- "D:/CodexWork/mic-re/meta/primer_table.csv"
cat("file exists:", file.exists(p), "\n")
cat("cutadapt exists:", file.exists("D:/bioenv/envs/riz16s/Scripts/cutadapt.exe"), "\n")
pt <- tryCatch(read.csv(p, stringsAsFactors=FALSE, check.names=FALSE, fileEncoding="UTF-8-BOM"),
               error=function(e) { cat("BOM read failed:", conditionMessage(e), "\n"); read.csv(p, stringsAsFactors=FALSE, check.names=FALSE) })
names(pt) <- sub("^\ufeff", "", names(pt))
cat("names:", paste(names(pt), collapse=" | "), "\n")
cat("nrow:", nrow(pt), "\n")
cat("projects:", paste(pt$project, collapse=", "), "\n")
proj <- "PRJNA1331937"
r <- pt[pt$project == proj, , drop=FALSE]
cat("matched rows for", proj, ":", nrow(r), "\n")
if (nrow(r)) { cat("primer_f =", r$primer_f[1], "\n"); cat("nzchar =", nzchar(r$primer_f[1]), "\n"); cat("is.na =", is.na(r$primer_f[1]), "\n") }
