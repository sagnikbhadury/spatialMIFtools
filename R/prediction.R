.patient_fold_features <- function(grids, phenotypes, neighborhood_fit,
                                   static_features = NULL) {
  assigned <- assign_neighborhoods(neighborhood_fit, grids)
  density <- paste0("density__", phenotypes)
  columns <- if (all(density %in% names(assigned))) density else phenotypes
  comp <- stats::aggregate(assigned[columns],
                           assigned[c("patient", "group")], mean, na.rm = TRUE)
  names(comp)[-(1:2)] <- paste0("composition__", phenotypes)
  nh <- as.data.frame.matrix(prop.table(table(assigned$patient,
                                              assigned$neighborhood), 1L))
  nh$patient <- rownames(nh)
  rownames(nh) <- NULL
  names(nh)[names(nh) != "patient"] <- paste0("neighborhood__",
                                                names(nh)[names(nh) != "patient"])
  out <- merge(comp, nh, by = "patient", all.x = TRUE)
  if (!is.null(static_features)) {
    .check_columns(static_features, "patient", "static_features")
    extra <- setdiff(names(static_features), c("group"))
    out <- merge(out, static_features[extra], by = "patient", all.x = TRUE)
    added <- setdiff(extra, "patient")
    names(out)[match(added, names(out))] <- paste0("static__", added)
  }
  out[is.na(out)] <- 0
  out
}

#' Fully nested leave-one-patient-out prediction
#'
#' Density scaling, K selection, K-means centers, held-out assignments, feature
#' scaling, and ridge-logistic fitting are repeated inside every training fold.
#'
#' @param grids Spatial grid table.
#' @param phenotypes Phenotypes used for composition and neighborhoods.
#' @param group_levels Ordered c(negative, positive) pathology groups.
#' @param feature_set composition, neighborhood, combined, or static.
#' @param static_features Optional one-row-per-patient features, such as network summaries.
#' @param k_range Candidate neighborhood numbers.
#' @param lambda Ridge-logistic penalty.
#' @param seed Random seed.
#' @return Out-of-fold patient scores with AUC as an attribute.
#' @export
nested_lopo_prediction <- function(grids, phenotypes, group_levels = NULL,
                                   feature_set = c("combined", "composition",
                                                   "neighborhood", "static"),
                                   static_features = NULL, k_range = 3:7,
                                   lambda = 4, seed = 20260803L) {
  feature_set <- match.arg(feature_set)
  .check_columns(grids, c("patient", "group"), "grids")
  patients <- unique(as.character(grids$patient))
  patient_group <- tapply(as.character(grids$group), grids$patient,
                          function(z) unique(z)[1L])
  groups <- unique(unname(patient_group))
  group_levels <- group_levels %||% groups[seq_len(min(2L, length(groups)))]
  if (length(group_levels) != 2L) stop("Exactly two group levels are required.", call. = FALSE)
  if (feature_set == "static" && is.null(static_features)) {
    stop("feature_set = 'static' requires static_features.", call. = FALSE)
  }
  rows <- vector("list", length(patients))
  for (fold in seq_along(patients)) {
    held <- patients[fold]
    train_grid <- grids[grids$patient != held, ]
    test_grid <- grids[grids$patient == held, ]
    fit <- fit_spatial_neighborhoods(train_grid, phenotypes, k_range = k_range,
                                     nstart = 10L, max_silhouette_bins = 2500L,
                                     seed = seed + fold)
    static_train <- if (is.null(static_features)) NULL else
      static_features[static_features$patient != held, , drop = FALSE]
    static_test <- if (is.null(static_features)) NULL else
      static_features[static_features$patient == held, , drop = FALSE]
    train <- .patient_fold_features(train_grid, phenotypes, fit, static_train)
    test <- .patient_fold_features(test_grid, phenotypes, fit, static_test)
    prefixes <- switch(feature_set,
                       composition = "composition__",
                       neighborhood = "neighborhood__",
                       static = "static__",
                       combined = c("composition__", "neighborhood__", "static__"))
    columns <- names(train)[vapply(names(train), function(z)
      any(startsWith(z, prefixes)), logical(1))]
    columns <- columns[vapply(train[columns], function(z) stats::sd(z) > 0, logical(1))]
    if (!length(columns)) stop("No varying predictors in fold ", fold, ".", call. = FALSE)
    x_train <- as.matrix(train[columns])
    x_test <- as.matrix(test[columns])
    scaler <- .safe_scale(x_train)
    x_train <- scaler$x[, scaler$active, drop = FALSE]
    x_test <- .apply_scaler(x_test, scaler)[, scaler$active, drop = FALSE]
    y_train <- as.integer(train$group == group_levels[2L])
    model <- .ridge_logistic_fit(x_train, y_train, lambda)
    probability <- .ridge_logistic_predict(model, x_test)
    rows[[fold]] <- data.frame(patient = held,
                               group = as.character(test$group[1L]),
                               observed = as.integer(test$group[1L] == group_levels[2L]),
                               probability = probability[1L], selected_k = fit$k,
                               n_training_patients = nrow(train),
                               n_predictors = ncol(x_train), stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  attr(out, "auc") <- .auc(out$observed, out$probability)
  attr(out, "group_levels") <- group_levels
  out
}
