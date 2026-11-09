# Training labels, gene-selection exclusions and input checks of the 20-gene classifier.
classifier_secondary_family <- "SF_INTESTINAL_GSE238274"
# Gene selection excludes only the four markers of the four-marker EC classifier. The 11
# primary-study-annotation genes (CCK, CHGA, GHRL, GIP, MLN, NEUROG3, NPW, NTS, PYY, SST, TAC1) are
# eligible everywhere the exclusion list is used: gene ranking, the neighborhood graph of the
# local-enrichment check, the matched-random draws and the recurrence null.
classifier_selection_exclusions <- c("FEV", "LMX1A", "SLC18A1", "TPH1")
# Genes are ranked by a study-weighted limma model (cp_training_gene_summary), the redundancy filter
# is off, and a panel classifies cells the same way as the four-marker classifier it is compared with:
# score threshold, gene detection and the local-enrichment check, with no threshold relaxation.
classifier_training_label_id <- "PRIMARY_LABEL_FIRST_THEN_FOUR_MARKER"
# Training labels: in the studies whose own deposited object names cells EC or a terminal non-EC
# enteroendocrine type (Egozi 2021, Elmentaite 2021, Hickey 2023), only those primary study labels
# are used and every other cell of the study abstains. The other studies (those without terminal
# primary study labels, and Elmentaite 2020, whose EC labels are a later marker-informed
# reannotation rather than the authors' own) take the four-marker EC classifier's label
# (METHOD_DEFINED_CONSERVATIVE in derive_comparison_labels).
classifier_primary_first_labels <- function(metadata) {
  require_columns(metadata, c("author_label_evidence_layer", "author_group_token"), "training metadata")
  rule <- derive_comparison_labels(metadata, "METHOD_DEFINED_CONSERVATIVE")
  flags <- population_flags(metadata)
  own <- metadata$author_label_evidence_layer %in% "ORIGINAL_AUTHOR_OBJECT"
  positive <- own & flags$exact_author_ec
  negative <- own & flags$exact_author_terminal_non_ec
  require_true(!any(positive & negative), "primary-label EC and non-EC masks overlap")
  label <- rule$label; reason <- paste0("AUTHOR_SILENT__", rule$reason)
  # In a study whose authors labeled terminal types, a cell they left broad or unresolved is not
  # given the four-marker classifier's label: it abstains, so that study is compared on primary
  # study labels alone.
  author_study <- metadata$family_id %in% unique(metadata$family_id[positive | negative])
  label[author_study] <- "ABSTAIN"; reason[author_study] <- "AUTHOR_STUDY_CELL_WITHOUT_TERMINAL_AUTHOR_LABEL"
  label[positive] <- "POSITIVE"; reason[positive] <- "ORIGINAL_AUTHOR_EC"
  label[negative] <- "NEGATIVE"; reason[negative] <- "ORIGINAL_AUTHOR_TERMINAL_NON_EC_EEC"
  data.frame(project_cell_id = metadata$project_cell_id, analysis_id = classifier_training_label_id,
             label = label, reason = reason, continuous_outcome = NA_real_, stringsAsFactors = FALSE)
}
# The training labels, after checking that every leave-one-study-out training and test split, and
# every inner split of the remaining studies, holds both classes.
classifier_preflight_labels <- function(metadata) {
  labels <- classifier_primary_first_labels(metadata)
  families <- sort(unique(metadata$family_id))
  assert_two <- function(keep, description) {
    require_true(all(c("POSITIVE", "NEGATIVE") %in% labels$label[keep]),
                 paste0(description, " lacks both classes"))
  }
  for (outer in families) {
    assert_two(metadata$family_id == outer, paste("outer test", outer))
    assert_two(metadata$family_id != outer, paste("outer training", outer))
    for (inner in setdiff(families, outer)) {
      assert_two(metadata$family_id == inner, paste("inner test", outer, inner))
      assert_two(!metadata$family_id %in% c(outer, inner), paste("inner training", outer, inner))
    }
  }
  assert_population_hierarchy(metadata)
  labels
}
classifier_assert_training_data <- function(x) {
  stopifnot(identical(x$part2_role, "TRAINING"),
            identical(x$part2_configuration, "egozi_included"),
            !any(x$metadata$family_id == classifier_secondary_family),
            any(grepl("EGOZI", x$metadata$family_id)),
            identical(colnames(x$normalized), x$metadata$project_cell_id),
            identical(colnames(x$counts), x$metadata$project_cell_id))
  invisible(TRUE)
}
classifier_assert_secondary_data <- function(x) {
  stopifnot(identical(x$part2_role, "SECONDARY_VALIDATION"),
            all(x$metadata$family_id == classifier_secondary_family),
            identical(colnames(x$normalized), x$metadata$project_cell_id))
  invisible(TRUE)
}
