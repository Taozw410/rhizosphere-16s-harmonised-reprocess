options(repos = c(CRAN = "https://cloud.r-project.org"), download.file.method = "libcurl", timeout = 120)
cat("BiocManager installed:", requireNamespace("BiocManager", quietly = TRUE), "\n")
u <- "https://cloud.r-project.org/src/contrib/PACKAGES.gz"
f <- tempfile()
r <- try(download.file(u, f, quiet = TRUE), silent = TRUE)
cat("CRAN reachable:", !inherits(r, "try-error"), " size:", if (file.exists(f)) file.info(f)$size else NA, "\n")
u2 <- "https://bioconductor.org/packages/release/bioc/src/contrib/PACKAGES.gz"
f2 <- tempfile()
r2 <- try(download.file(u2, f2, quiet = TRUE), silent = TRUE)
cat("Bioconductor reachable:", !inherits(r2, "try-error"), " size:", if (file.exists(f2)) file.info(f2)$size else NA, "\n")
cat("R version:", R.version.string, "\n")
