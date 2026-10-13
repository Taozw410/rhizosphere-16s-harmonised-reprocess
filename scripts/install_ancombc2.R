z <- "D:/CodexWork/mic-re/out_analysis/pkgzips/ANCOMBC_2.10.1.zip"
r <- try(install.packages(z, repos = NULL, type = "binary"), silent = TRUE)
cat("ok:", requireNamespace("ANCOMBC", quietly = TRUE), "\n")
if (!requireNamespace("ANCOMBC", quietly = TRUE)) {
  cat("error:", substr(as.character(r), 1, 200), "\n")
  ap <- available.packages(type = "binary")
  d <- tools::package_dependencies("ANCOMBC", db = ap, recursive = TRUE)[[1]]
  inst <- rownames(installed.packages())
  cat("still missing:", paste(setdiff(d, inst), collapse = ", "), "\n")
}
