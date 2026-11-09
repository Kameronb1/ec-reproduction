# Reads an .h5ad file (CSR X, optional CSR layers, obs/var with plain, categorical and nullable
# columns) into a SingleCellExperiment with rhdf5 only, as zellkonverter::readH5AD(use_hdf5 = FALSE)
# would. It stops on any other encoding. Sourced by heldout.R and four_marker_depth.R.

suppressPackageStartupMessages({
  library(rhdf5)
  library(Matrix)
  library(SingleCellExperiment)
  library(S4Vectors)
})

h5ad_fail <- function(...) stop(paste0("READ_H5AD_RHDF5_FAIL: ", ...), call. = FALSE)

# anndata writes strings as variable-length UTF-8; rhdf5 returns them as character already.
.h5ad_chr <- function(x) {
  if (is.raw(x)) h5ad_fail("raw string payload is not supported")
  as.character(x)
}

# One obs/var column. Three encodings appear in these objects and nothing else is accepted.
.h5ad_column <- function(value, name) {
  if (!is.list(value)) {
    return(if (is.character(value)) .h5ad_chr(value) else as.vector(value))
  }
  members <- names(value)
  if (all(c("categories", "codes") %in% members)) {
    categories <- .h5ad_chr(as.vector(value$categories))
    codes <- as.integer(as.vector(value$codes))
    if (any(codes < -1L) || any(codes >= length(categories))) {
      h5ad_fail("column ", name, " has a categorical code outside its categories")
    }
    out <- rep(NA_character_, length(codes))
    present <- codes >= 0L
    out[present] <- categories[codes[present] + 1L]
    # zellkonverter surfaces pandas categoricals as factors with the stored category order.
    return(factor(out, levels = categories))
  }
  if ("values" %in% members) {
    values <- as.vector(value$values)
    if (is.character(values)) values <- .h5ad_chr(values)
    if ("mask" %in% members) {
      mask <- as.logical(as.vector(value$mask))
      if (length(mask) == length(values)) values[mask] <- NA
    }
    return(values)
  }
  h5ad_fail("column ", name, " uses an unsupported encoding: ", paste(members, collapse = ", "))
}

# An obs or var group as a data frame, with the anndata _index kept both as the row names and as
# a column, as zellkonverter keeps it.
.h5ad_frame <- function(path, group, expected_rows) {
  raw <- rhdf5::h5read(path, group)
  attrs <- rhdf5::h5readAttributes(path, group)
  index_name <- attrs[["_index"]]
  if (is.null(index_name)) h5ad_fail(group, " has no _index attribute")
  index_name <- .h5ad_chr(index_name)[[1L]]
  if (!(index_name %in% names(raw))) h5ad_fail(group, " index column is absent: ", index_name)

  column_order <- attrs[["column-order"]]
  ordered <- if (!is.null(column_order)) .h5ad_chr(as.vector(column_order)) else character(0)
  ordered <- ordered[ordered %in% names(raw)]
  remaining <- setdiff(names(raw), c(ordered, index_name))
  keys <- c(index_name, ordered, remaining)

  columns <- lapply(keys, function(key) .h5ad_column(raw[[key]], paste0(group, "/", key)))
  names(columns) <- keys
  lengths_seen <- vapply(columns, length, integer(1))
  if (any(lengths_seen != expected_rows)) {
    h5ad_fail(group, " column lengths disagree with the matrix: ",
              paste(unique(lengths_seen), collapse = ", "), " against ", expected_rows)
  }
  frame <- as.data.frame(columns, stringsAsFactors = FALSE, check.names = FALSE)
  rownames(frame) <- as.character(columns[[index_name]])
  frame
}

# anndata stores X as CSR over CELLS. An R dgCMatrix is CSC, so a dgCMatrix whose Dim is
# (features, cells) has exactly the same three arrays: the column pointer is the CSR row pointer
# and the row index is the CSR column index. The transpose to features x cells is therefore free
# and exact - no data is moved, reordered or recomputed.
.h5ad_matrix <- function(path, group, n_obs, n_vars, obs_names, var_names) {
  attrs <- rhdf5::h5readAttributes(path, group)
  encoding <- attrs[["encoding-type"]]
  if (is.null(encoding) || !identical(.h5ad_chr(encoding)[[1L]], "csr_matrix")) {
    h5ad_fail(group, " is not a csr_matrix; this reader does not guess at other encodings")
  }
  shape <- as.integer(as.vector(attrs[["shape"]]))
  if (!identical(shape, c(as.integer(n_obs), as.integer(n_vars)))) {
    h5ad_fail(group, " shape is ", paste(shape, collapse = "x"),
              ", expected ", n_obs, "x", n_vars)
  }
  indptr <- as.vector(rhdf5::h5read(path, paste0(group, "/indptr")))
  indices <- as.vector(rhdf5::h5read(path, paste0(group, "/indices")))
  data <- as.vector(rhdf5::h5read(path, paste0(group, "/data")))
  if (length(indptr) != n_obs + 1L) h5ad_fail(group, " indptr length is wrong")
  if (length(indices) != length(data)) h5ad_fail(group, " indices and data lengths differ")
  if (indptr[[1L]] != 0) h5ad_fail(group, " indptr does not start at 0")
  if (indptr[[length(indptr)]] != length(data)) h5ad_fail(group, " indptr does not end at nnz")
  if (length(indices) && (min(indices) < 0 || max(indices) >= n_vars)) {
    h5ad_fail(group, " has a column index outside the feature range")
  }
  matrix <- new(
    "dgCMatrix",
    i = as.integer(indices),
    p = as.integer(indptr),
    x = as.numeric(data),
    Dim = c(as.integer(n_vars), as.integer(n_obs)),
    Dimnames = list(as.character(var_names), as.character(obs_names))
  )
  validObject(matrix)
  matrix
}

#' Read an .h5ad into a SingleCellExperiment the way readH5AD(use_hdf5 = FALSE) would.
read_h5ad_as_sce <- function(path) {
  if (!file.exists(path)) h5ad_fail("file is absent: ", path)
  contents <- rhdf5::h5ls(path, recursive = 1L)
  top <- contents$name[contents$group == "/"]
  for (required in c("X", "obs", "var")) {
    if (!(required %in% top)) h5ad_fail("required group is absent: ", required)
  }
  x_attrs <- rhdf5::h5readAttributes(path, "X")
  shape <- as.integer(as.vector(x_attrs[["shape"]]))
  if (length(shape) != 2L) h5ad_fail("X has no two-dimensional shape attribute")
  n_obs <- shape[[1L]]
  n_vars <- shape[[2L]]

  var_frame <- .h5ad_frame(path, "var", n_vars)
  obs_frame <- .h5ad_frame(path, "obs", n_obs)
  var_names <- rownames(var_frame)
  obs_names <- rownames(obs_frame)

  assays <- list(X = .h5ad_matrix(path, "X", n_obs, n_vars, obs_names, var_names))
  if ("layers" %in% top) {
    # recursive = 1 stops before the children of /layers, so the layer names need depth 2.
    listing <- rhdf5::h5ls(path, recursive = 2L)
    layer_names <- unique(listing$name[listing$group == "/layers"])
    if (!length(layer_names)) h5ad_fail("the layers group is present but lists no layer")
    for (layer in layer_names) {
      assays[[layer]] <- .h5ad_matrix(path, paste0("layers/", layer),
                                      n_obs, n_vars, obs_names, var_names)
    }
  }
  rhdf5::h5closeAll()

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = assays,
    rowData = S4Vectors::DataFrame(var_frame, check.names = FALSE),
    colData = S4Vectors::DataFrame(obs_frame, check.names = FALSE)
  )
  rownames(sce) <- var_names
  colnames(sce) <- obs_names
  sce
}
