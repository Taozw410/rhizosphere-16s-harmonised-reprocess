z <- "D:/CodexWork/mic-re/out_analysis/pkgzips/ANCOMBC_2.10.1.zip"
cat("zip exists:", file.exists(z), " size KB:", round(file.info(z)$size/1024,1), "\n")
r <- try(install.packages(z, repos = NULL, type = "binary"), silent = FALSE)
cat("class:", class(r), "\n")
cat("ANCOMBC available:", requireNamespace("ANCOMBC", quietly = TRUE), "\n")
if (!requireNamespace("ANCOMBC", quietly = TRUE)) {
  # which of its hard deps are missing?
  ap <- available.packages(type="binary")
  d <- tools::package_dependencies("ANCOMBC", db = ap, recursive = FALSE)[[1]]
  inst <- rownames(installed.packages())
  cat("missing hard deps:", paste(setdiff(d, inst), collapse=", "), "\n")
}
