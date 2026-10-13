options(repos = c(CRAN="https://cloud.r-project.org",
                  BioCsoft="https://bioconductor.org/packages/3.21/bioc",
                  BioCann="https://bioconductor.org/packages/3.21/data/annotation",
                  BioCexp="https://bioconductor.org/packages/3.21/data/experiment"), timeout = 3000)
ap <- available.packages(type = "binary")
d <- tools::package_dependencies("ANCOMBC", db = ap, recursive = TRUE)[[1]]
inst <- rownames(installed.packages())
miss <- setdiff(d, inst)
cat("missing deps:", length(miss), "\n")
cat(paste(miss, collapse = ", "), "\n")
rows <- data.frame()
for (p in miss) if (p %in% rownames(ap))
  rows <- rbind(rows, data.frame(pkg = p, ver = ap[p, "Version"], url = ap[p, "Repository"], stringsAsFactors = FALSE))
write.csv(rows, "D:/CodexWork/mic-re/out_analysis/install_plan2.csv", row.names = FALSE)
cat("planned:", nrow(rows), "\n")
print(rows)
