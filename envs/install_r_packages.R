# Install the R packages the analysis scripts use, at the versions the published tables were
# checked with (R 4.6.1, Bioconductor 3.23). Run once:  Rscript envs/install_r_packages.R
# Then run envs/check_versions.R; it lists any package whose installed version differs.
stopifnot(getRversion() >= "4.6.0")
pkgs <- read.delim(file.path("envs", "r_packages.tsv"), stringsAsFactors = FALSE)
options(repos = c(CRAN = "https://cloud.r-project.org"))
for (p in c("remotes", "BiocManager")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
BiocManager::install(version = "3.23", ask = FALSE, update = FALSE)
cran <- pkgs[pkgs$source == "CRAN", ]
for (i in seq_len(nrow(cran))) {
  have <- tryCatch(as.character(packageVersion(cran$package[i])), error = function(e) "")
  if (have != cran$version[i]) remotes::install_version(cran$package[i], cran$version[i], upgrade = "never")
}
bioc <- pkgs$package[pkgs$source != "CRAN"]
BiocManager::install(bioc, ask = FALSE, update = FALSE)
