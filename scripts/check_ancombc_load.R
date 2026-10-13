lib <- "D:/bioenv/envs/riz16s/Lib/R/library"
d <- file.path(lib, "ANCOMBC")
cat("dir exists:", dir.exists(d), "\n")
if (dir.exists(d)) cat("files:", paste(head(list.files(d), 8), collapse=", "), "\n")
cat("DESCRIPTION exists:", file.exists(file.path(d, "DESCRIPTION")), "\n")
cat("--- loading ---\n")
r <- try(library(ANCOMBC), silent = TRUE)
if (inherits(r, "try-error")) cat("LOAD ERROR:", substr(as.character(attr(r, "condition")$message), 1, 300), "\n")
cat("loaded:", "ANCOMBC" %in% loadedNamespaces(), "\n")
