# List R packages whose installed version differs from envs/r_packages.tsv.
# A difference does not stop the workflow; the table checks decide whether results match.
pkgs <- read.delim(file.path("envs", "r_packages.tsv"), stringsAsFactors = FALSE)
have <- vapply(pkgs$package, function(p) tryCatch(as.character(packageVersion(p)), error = function(e) "missing"), "")
diff <- pkgs$version != have
cat("R", as.character(getRversion()), "(tables checked with 4.6.1)\n")
if (any(diff)) print(data.frame(package = pkgs$package, wanted = pkgs$version, installed = have)[diff, ], row.names = FALSE) else cat("all", nrow(pkgs), "packages match\n")
