#' Ridge partial-correlation matrix
#' @param x Numeric observations by phenotype matrix.
#' @param ridge Positive diagonal ridge penalty.
#' @return Partial-correlation matrix; constant variables are non-estimable.
#' @export
ridge_partial_correlation <- function(x, ridge = 0.15) {
  x <- as.matrix(x)
  if (!is.numeric(x) || nrow(x) < 2L) stop("x must be a numeric matrix.", call. = FALSE)
  if (!is.finite(ridge) || ridge < 0) stop("ridge must be nonnegative.", call. = FALSE)
  s <- apply(x, 2L, stats::sd, na.rm = TRUE)
  active <- is.finite(s) & s > 0
  out <- matrix(NA_real_, ncol(x), ncol(x), dimnames = list(colnames(x), colnames(x)))
  diag(out) <- 1
  if (sum(active) < 2L) return(out)
  z <- scale(x[, active, drop = FALSE])
  covariance <- crossprod(z) / max(nrow(z) - 1L, 1L) + ridge * diag(sum(active))
  precision <- solve(covariance)
  pcor <- -precision / sqrt(outer(diag(precision), diag(precision)))
  diag(pcor) <- 1
  out[active, active] <- pcor
  out
}

#' Residualize phenotype densities on local total cellularity
#' @param grids Grid data.
#' @param phenotypes Phenotype names or density-column names.
#' @param total_column Total cell-count column.
#' @return Numeric residual matrix.
#' @export
residualize_cellularity <- function(grids, phenotypes,
                                    total_column = "n_cells") {
  density <- paste0("density__", phenotypes)
  columns <- if (all(density %in% names(grids))) density else phenotypes
  .check_columns(grids, c(columns, total_column), "grids")
  total <- log1p(grids[[total_column]])
  out <- vapply(columns, function(column) {
    stats::residuals(stats::lm(grids[[column]] ~ total))
  }, numeric(nrow(grids)))
  colnames(out) <- phenotypes
  out
}

#' Estimate image-level ridge conditional networks
#' @param grids Spatial grids.
#' @param phenotypes Phenotype names.
#' @param ridge Ridge penalty.
#' @param minimum_bins Minimum bins per image or image-zone.
#' @param adjust_cellularity Residualize on local total count first.
#' @param by_zone Estimate networks separately within tumor_zone.
#' @return Long-form image-level edge table.
#' @export
estimate_ridge_networks <- function(grids, phenotypes, ridge = 0.15,
                                    minimum_bins = 12L,
                                    adjust_cellularity = FALSE,
                                    by_zone = FALSE) {
  required <- c("patient", "image", "group")
  if (by_zone) required <- c(required, "tumor_zone")
  .check_columns(grids, required, "grids")
  minimum_bins <- .assert_positive_integer(minimum_bins, "minimum_bins")
  split_key <- interaction(grids$patient, grids$image,
                           if (by_zone) grids$tumor_zone else "overall", drop = TRUE)
  pieces <- split(grids, split_key)
  pairs <- .pair_names(phenotypes)
  rows <- lapply(pieces, function(d) {
    if (nrow(d) < minimum_bins) return(NULL)
    density <- paste0("density__", phenotypes)
    columns <- if (all(density %in% names(d))) density else phenotypes
    .check_columns(d, columns, "grids")
    x <- if (adjust_cellularity) residualize_cellularity(d, phenotypes) else {
      ans <- as.matrix(d[columns]); colnames(ans) <- phenotypes; ans
    }
    pcor <- ridge_partial_correlation(x, ridge)
    out <- pairs
    out$partial_correlation <- mapply(function(a, b) pcor[a, b], out$source, out$target)
    out$patient <- as.character(d$patient[1L])
    out$image <- as.character(d$image[1L])
    out$group <- as.character(d$group[1L])
    out$zone <- if (by_zone) as.character(d$tumor_zone[1L]) else "overall"
    out$n_bins <- nrow(d)
    out$cellularity_adjusted <- adjust_cellularity
    out
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Compare patient-averaged conditional edges between two groups
#' @param networks Output from [estimate_ridge_networks()].
#' @param group_levels Ordered c(reference, comparison) group names.
#' @param p_adjust Multiplicity procedure passed to [stats::p.adjust()].
#' @param bootstrap_replicates Patient bootstrap replicates for percentile intervals.
#' @param seed Random seed.
#' @return One row per edge with comparison-minus-reference estimates.
#' @export
patient_edge_contrasts <- function(networks, group_levels = NULL,
                                   p_adjust = "BH", bootstrap_replicates = 0L,
                                   seed = 20260803L) {
  .check_columns(networks, c("patient", "group", "edge", "partial_correlation"),
                 "networks")
  patient <- stats::aggregate(partial_correlation ~ patient + group + edge,
                              networks, mean, na.rm = TRUE)
  groups <- unique(as.character(patient$group))
  group_levels <- group_levels %||% groups[seq_len(min(2L, length(groups)))]
  if (length(group_levels) != 2L || !all(group_levels %in% groups)) {
    stop("group_levels must name exactly two observed groups.", call. = FALSE)
  }
  edges <- unique(patient$edge)
  rows <- lapply(edges, function(edge) {
    d <- patient[patient$edge == edge, ]
    a <- d$partial_correlation[d$group == group_levels[1L]]
    b <- d$partial_correlation[d$group == group_levels[2L]]
    test <- stats::t.test(b, a, var.equal = FALSE)
    ci <- c(NA_real_, NA_real_)
    if (bootstrap_replicates > 0L) {
      set.seed(seed + match(edge, edges))
      boots <- replicate(bootstrap_replicates,
                         mean(sample(b, length(b), TRUE)) - mean(sample(a, length(a), TRUE)))
      ci <- stats::quantile(boots, c(.025, .975), na.rm = TRUE, names = FALSE)
    }
    data.frame(edge = edge, reference = group_levels[1L], comparison = group_levels[2L],
               mean_reference = mean(a), mean_comparison = mean(b),
               estimate = mean(b) - mean(a), statistic = unname(test$statistic),
               p_value = test$p.value, conf_low = ci[1L], conf_high = ci[2L])
  })
  out <- do.call(rbind, rows)
  out$q_value <- stats::p.adjust(out$p_value, method = p_adjust)
  out[order(out$q_value, out$p_value), ]
}
