#' Four-nearest-neighbor spatial autocorrelation score
#'
#' The statistic is proportional to Moran's I on a four-neighbor graph but is
#' not normalized by total spatial weight.
#'
#' @param values Numeric vector at occupied bins.
#' @param coords Two-column coordinate matrix.
#' @param permutations Number of within-image randomizations.
#' @param seed Random seed.
#' @return A list with score and permutation p_value.
#' @export
four_nn_autocorrelation <- function(values, coords, permutations = 199L,
                                    seed = NULL) {
  values <- as.numeric(values)
  coords <- as.matrix(coords)
  if (nrow(coords) != length(values) || ncol(coords) != 2L) {
    stop("coords must have two columns and one row per value.", call. = FALSE)
  }
  ok <- is.finite(values) & stats::complete.cases(coords)
  values <- values[ok]
  coords <- coords[ok, , drop = FALSE]
  if (length(values) < 5L || stats::sd(values) == 0) {
    return(list(score = NA_real_, p_value = NA_real_))
  }
  d <- as.matrix(stats::dist(coords))
  diag(d) <- Inf
  neighbors <- t(apply(d, 1L, order))[, 1:4, drop = FALSE]
  score_one <- function(v) {
    z <- v - mean(v)
    sum(vapply(seq_along(z), function(i) sum(z[i] * z[neighbors[i, ]]),
               numeric(1))) / sum(z^2)
  }
  observed <- score_one(values)
  permutations <- as.integer(permutations)
  if (permutations < 1L) return(list(score = observed, p_value = NA_real_))
  if (!is.null(seed)) set.seed(seed)
  null <- replicate(permutations, score_one(sample(values)))
  p <- (1 + sum(abs(null) >= abs(observed))) / (permutations + 1)
  list(score = observed, p_value = p)
}

.silhouette_mean <- function(x, clusters, max_n = 6000L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (nrow(x) > max_n) {
    keep <- sample.int(nrow(x), max_n)
    x <- x[keep, , drop = FALSE]
    clusters <- clusters[keep]
  }
  d <- as.matrix(stats::dist(x))
  vals <- vapply(seq_len(nrow(x)), function(i) {
    own <- which(clusters == clusters[i] & seq_len(nrow(x)) != i)
    a <- if (length(own)) mean(d[i, own]) else 0
    others <- setdiff(unique(clusters), clusters[i])
    b <- min(vapply(others, function(k) mean(d[i, clusters == k]), numeric(1)))
    if (max(a, b) == 0) 0 else (b - a) / max(a, b)
  }, numeric(1))
  mean(vals)
}

#' Discover recurrent multivariate spatial neighborhoods
#'
#' @param grids Spatial grid table.
#' @param phenotypes Phenotype names; density__ columns are used when present.
#' @param k_range Candidate neighborhood numbers.
#' @param nstart K-means random initializations.
#' @param max_silhouette_bins Maximum bins used for silhouette evaluation.
#' @param seed Random seed.
#' @return A fitted spatialmif_neighborhoods object.
#' @export
fit_spatial_neighborhoods <- function(grids, phenotypes,
                                      k_range = 3:7, nstart = 50L,
                                      max_silhouette_bins = 6000L,
                                      seed = 20260803L) {
  density <- paste0("density__", phenotypes)
  columns <- if (all(density %in% names(grids))) density else phenotypes
  .check_columns(grids, columns, "grids")
  x <- as.matrix(grids[columns])
  scaler <- .safe_scale(x)
  active <- scaler$active
  if (sum(active) < 1L) stop("No varying phenotype density is available.", call. = FALSE)
  z <- scaler$x[, active, drop = FALSE]
  k_range <- sort(unique(as.integer(k_range)))
  k_range <- k_range[k_range >= 2L & k_range < nrow(z)]
  if (!length(k_range)) stop("k_range has no admissible value.", call. = FALSE)
  set.seed(seed)
  fits <- lapply(k_range, function(k) stats::kmeans(z, centers = k, nstart = nstart))
  scores <- vapply(seq_along(fits), function(i) {
    .silhouette_mean(z, fits[[i]]$cluster, max_silhouette_bins, seed + i)
  }, numeric(1))
  best <- which.max(scores)
  fit <- fits[[best]]
  centers <- matrix(0, nrow(fit$centers), length(columns),
                    dimnames = list(paste0("N", seq_len(nrow(fit$centers))), columns))
  centers[, active] <- fit$centers
  structure(list(k = k_range[best], centers = centers, scaler = scaler,
                 columns = columns, phenotypes = phenotypes,
                 scores = data.frame(k = k_range, silhouette = scores),
                 cluster = fit$cluster), class = "spatialmif_neighborhoods")
}

#' Assign bins to fitted spatial neighborhoods
#' @param object A fitted object from [fit_spatial_neighborhoods()].
#' @param grids Grid table to assign.
#' @param append Return the grid with a neighborhood column when TRUE.
#' @return Integer assignments or the augmented grid.
#' @export
assign_neighborhoods <- function(object, grids, append = TRUE) {
  if (!inherits(object, "spatialmif_neighborhoods")) {
    stop("object must come from fit_spatial_neighborhoods().", call. = FALSE)
  }
  .check_columns(grids, object$columns, "grids")
  z <- .apply_scaler(as.matrix(grids[object$columns]), object$scaler)
  dist_to_center <- vapply(seq_len(nrow(object$centers)), function(k) {
    rowSums((z - matrix(object$centers[k, ], nrow(z), ncol(z), byrow = TRUE))^2)
  }, numeric(nrow(z)))
  assignment <- max.col(-dist_to_center, ties.method = "first")
  if (!append) return(assignment)
  grids$neighborhood <- paste0("N", assignment)
  grids
}

#' Patient-level neighborhood prevalence
#' @param grids Assigned grid returned by [assign_neighborhoods()].
#' @return Long-form patient-level prevalence table.
#' @export
patient_neighborhood_prevalence <- function(grids) {
  .check_columns(grids, c("patient", "group", "neighborhood"), "grids")
  tab <- as.data.frame(table(patient = grids$patient, group = grids$group,
                             neighborhood = grids$neighborhood),
                       stringsAsFactors = FALSE)
  tab <- tab[tab$Freq > 0, ]
  totals <- stats::aggregate(Freq ~ patient + group, tab, sum)
  names(totals)[3L] <- "total"
  tab <- merge(tab, totals, by = c("patient", "group"), all.x = TRUE)
  tab$prevalence <- tab$Freq / tab$total
  tab[c("patient", "group", "neighborhood", "prevalence")]
}
