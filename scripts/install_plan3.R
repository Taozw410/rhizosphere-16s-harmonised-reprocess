options(repos = c(CRAN="https://cloud.r-project.org",
                  BioCsoft="https://bioconductor.org/packages/3.21/bioc",
                  BioCann="https://bioconductor.org/packages/3.21/data/annotation",
                  BioCexp="https://bioconductor.org/packages/3.21/data/experiment"), timeout=3600)
ap <- available.packages(type="binary")
cat("index rows:", nrow(ap), "\n")
deps <- tools::package_dependencies("microbiome", db=ap, recursive=TRUE)[[1]]
miss <- setdiff(c("microbiome", deps), rownames(installed.packages()))
rows <- data.frame()
for (p in miss) if (p %in% rownames(ap)) rows <- rbind(rows, data.frame(pkg=p, ver=ap[p,"Version"], url=ap[p,"Repository"], stringsAsFactors=FALSE))
write.csv(rows, "D:/CodexWork/mic-re/out_analysis/install_plan3.csv", row.names=FALSE)
cat("missing:", length(miss), " planned:", nrow(rows), "\n")
cat(paste(rows$pkg, collapse=", "), "\n")
