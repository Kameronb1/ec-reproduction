# Checks, the label comparisons and the gene-eligibility filter of the 20-gene classifier.
# Sourcing this file reads no data, fits no model and writes nothing.

abortf <- function(format, ...) {
  stop(sprintf(format, ...), call. = FALSE)
}

require_true <- function(value, message) {
  if (!isTRUE(value)) abortf("%s", message)
  invisible(TRUE)
}

require_columns <- function(frame, columns, label = "table") {
  missing <- setdiff(columns, names(frame))
  require_true(length(missing) == 0L, sprintf("%s missing columns: %s", label, paste(missing, collapse = ", ")))
  invisible(TRUE)
}

# The three classes of the four-marker EC classifier. Any other value fails before fitting.
CLASSIFIER_CLASSES <- c(high = "HIGH_CONFIDENCE_EC", low = "LOW_CONFIDENCE_EC",
                  non_ec = "INSUFFICIENT_EC_EVIDENCE")
population_flags <- function(metadata) {
  require_columns(metadata, c("presentation_token", "author_group_token", "ec_core"), "standardized metadata")
  presentation <- metadata$presentation_token
  require_true(!anyNA(presentation) && all(presentation %in% CLASSIFIER_CLASSES),
               paste0("unknown or missing four-marker class: ",
                      paste(setdiff(unique(presentation), CLASSIFIER_CLASSES), collapse = ", ")))
  author <- metadata$author_group_token
  data.frame(
    high_confidence_ec = presentation == CLASSIFIER_CLASSES[["high"]],
    low_confidence_ec = presentation == CLASSIFIER_CLASSES[["low"]],
    non_ec = presentation == CLASSIFIER_CLASSES[["non_ec"]],
    exact_author_ec = author %in% c("EXPLICIT_AUTHOR_EC", "EXPLICIT_EC", "EXACT_AUTHOR_EC", "AUTHOR_EC"),
    exact_author_terminal_non_ec = author %in% c("EXPLICIT_TERMINAL_NON_EC_EEC", "TERMINAL_NON_EC_EEC", "EXACT_AUTHOR_TERMINAL_NON_EC_EEC", "AUTHOR_TERMINAL_NON_EC_EEC"),
    stringsAsFactors = FALSE
  )
}

derive_comparison_labels <- function(metadata, analysis_id) {
  flags <- population_flags(metadata)
  label <- rep("ABSTAIN", nrow(metadata))
  reason <- rep("NOT_IN_COMPARISON", nrow(metadata))
  if (identical(analysis_id, "METHOD_DEFINED_CONSERVATIVE")) {
    positive <- flags$high_confidence_ec & !flags$exact_author_terminal_non_ec
    negative <- flags$non_ec
    require_true(!any(positive & negative), "method-defined comparison masks overlap")
    label[positive] <- "POSITIVE"
    reason[positive] <- "HIGH_CONFIDENCE_WITHOUT_AUTHOR_TERMINAL_CONFLICT"
    label[negative] <- "NEGATIVE"
    reason[negative] <- "RULE_NON_EC_CALL"
    reason[flags$low_confidence_ec] <- "LOW_CONFIDENCE_ABSTENTION"
    reason[flags$high_confidence_ec & flags$exact_author_terminal_non_ec] <- "HIGH_CONFIDENCE_AUTHOR_TERMINAL_CONFLICT"
  } else if (identical(analysis_id, "DUAL_POSITIVE_INTERSECTION")) {
    positive <- flags$high_confidence_ec & flags$exact_author_ec
    negative <- flags$exact_author_terminal_non_ec & flags$non_ec
    require_true(!any(positive & negative), "intersection comparison masks overlap")
    label[positive] <- "POSITIVE"
    reason[positive] <- "HIGH_CONFIDENCE_AND_EXACT_AUTHOR_EC"
    label[negative] <- "NEGATIVE"
    reason[negative] <- "EXACT_TERMINAL_NON_EC_AND_NOT_HIGH_CONFIDENCE"
  } else if (identical(analysis_id, "AUTHOR_LABEL_ONLY")) {
    require_true(!any(flags$exact_author_ec & flags$exact_author_terminal_non_ec), "primary-label comparison masks overlap")
    label[flags$exact_author_ec] <- "POSITIVE"
    reason[flags$exact_author_ec] <- "EXACT_AUTHOR_EC"
    label[flags$exact_author_terminal_non_ec] <- "NEGATIVE"
    reason[flags$exact_author_terminal_non_ec] <- "EXACT_AUTHOR_TERMINAL_NON_EC_EEC"
  } else if (identical(analysis_id, "CONTINUOUS_EC_CORE")) {
    return(data.frame(
      project_cell_id = metadata$project_cell_id,
      analysis_id = analysis_id,
      label = "CONTINUOUS",
      reason = "FIXED_EC_CORE_SCORE",
      continuous_outcome = metadata$ec_core,
      stringsAsFactors = FALSE
    ))
  } else {
    abortf("unknown comparison analysis_id: %s", analysis_id)
  }
  data.frame(
    project_cell_id = metadata$project_cell_id,
    analysis_id = analysis_id,
    label = label,
    reason = reason,
    continuous_outcome = NA_real_,
    stringsAsFactors = FALSE
  )
}

assert_population_hierarchy <- function(metadata) {
  primary <- derive_comparison_labels(metadata, "METHOD_DEFINED_CONSERVATIVE")
  intersection <- derive_comparison_labels(metadata, "DUAL_POSITIVE_INTERSECTION")
  author <- derive_comparison_labels(metadata, "AUTHOR_LABEL_ONLY")
  primary_positive <- primary$project_cell_id[primary$label == "POSITIVE"]
  intersection_positive <- intersection$project_cell_id[intersection$label == "POSITIVE"]
  author_positive <- author$project_cell_id[author$label == "POSITIVE"]
  require_true(all(intersection_positive %in% primary_positive), "intersection positives are not a subset of method-defined positives")
  require_true(all(intersection_positive %in% author_positive), "intersection positives are not a subset of primary-label positives")
  invisible(TRUE)
}

training_feature_mask <- function(log_expression, genes, metadata, train, detection_floor, minimum_training_families, excluded_genes) {
  require_true(inherits(log_expression, "Matrix") || is.matrix(log_expression), "log expression must be a matrix-like genes-by-cells object")
  require_true(nrow(log_expression) == length(genes) && ncol(log_expression) == nrow(metadata), "expression dimensions do not match genes and metadata")
  require_true(is.numeric(detection_floor) && length(detection_floor) == 1L && detection_floor > 0 && detection_floor < 1, "detection_floor must be in (0,1)")
  families <- sort(unique(metadata$family_id[train]))
  pass_count <- integer(length(genes))
  for (family in families) {
    cells <- train & metadata$family_id == family
    pass_count <- pass_count + (Matrix::rowMeans(log_expression[, cells, drop = FALSE] > 0) >= detection_floor)
  }
  eligible <- pass_count >= as.integer(minimum_training_families)
  eligible <- eligible & !(genes %in% excluded_genes)
  eligible <- eligible & !grepl("^MT-", genes, ignore.case = TRUE)
  eligible <- eligible & !grepl("^RP[SL][0-9]", genes, ignore.case = TRUE)
  eligible
}
