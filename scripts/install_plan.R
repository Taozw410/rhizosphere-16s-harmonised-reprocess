repos <- c(CRAN = "https://cloud.r-project.org",
           BioCsoft = "https://bioconductor.org/packages/3.21/bioc",
           BioCann  = "https://bioconductor.org/packages/3.21/data/annotation",
           BioCexp  = "https://bioconductor.org/packages/3.21/data/experiment")
options(repos = repos, download.file.method = "libcurl", timeout = 3600)
want <- c("ALDEx2", "ANCOMBC")
ap <- try(available.packages(type = "binary"), silent = TRUE)
if (inherits(ap, "try-error")) { cat("cannot fetch binary index\n"); quit(status = 1) }
cat("binary index rows:", nrow(ap), "\n")
dep <- try(tools::package_dependencies(want, db = ap, recursive = TRUE), silent = TRUE)
need <- unique(c(want, unlist(dep)))
inst <- rownames(installed.packages())
missing <- setdiff(need, inst)
cat("total in closure:", length(need), " already installed:", length(intersect(need, inst)), " missing:", length(missing), "\n")
rows <- data.frame()
for (p in missing) {
  if (!p %in% rownames(ap)) { cat("  NOT in binary index:", p, "\n"); next }
  rows <- rbind(rows, data.frame(pkg = p, ver = ap[p, "Version"], url = ap[p, "Repository"], stringsAsFactors = FALSE))
}
out <- "D:/CodexWork/mic-re/out_analysis/install_plan.csv"
write.csv(rows, out, row.names = FALSE)
cat("planned downloads:", nrow(rows), "->", out, "\n")
print(rows)
