.cross_pair_count <- function(distance, a, b, radius) {
  pair <- which(upper.tri(distance) & distance <= radius, arr.ind = TRUE)
  if (!nrow(pair)) return(0L)
  sum((a[pair[, 1L]] & b[pair[, 2L]]) |
        (b[pair[, 1L]] & a[pair[, 2L]]))
}

#' Multiradius random-label cross-pair enrichment
#'
#' @param cells Physical-cell table.
#' @param pairs Two-column matrix/data frame or list of phenotype pairs.
#' @param radii Distances after separately normalizing x and y by image span.
#' @param permutations Random relabelings within image.
#' @param seed Random seed.
#' @return Image-level observed, expected, enrichment, and random-label p-values.
#' @export
cross_pair_enrichment <- function(cells, pairs, radii = c(.025, .05, .10),
                                  permutations = 99L, seed = 20260803L) {
  .check_columns(cells, .required_cell_columns, "cells")
  if (is.list(pairs) && !is.data.frame(pairs)) pairs <- do.call(rbind, pairs)
  pairs <- as.matrix(pairs)
  if (ncol(pairs) != 2L) stop("pairs must have two phenotype columns.", call. = FALSE)
  .check_columns(cells, unique(as.vector(pairs)), "cells")
  permutations <- .assert_positive_integer(permutations, "permutations")
  set.seed(seed)
  pieces <- .image_split(cells)
  rows <- list()
  k <- 0L
  for (d in pieces) {
    xy <- as.matrix(d[c("x", "y")])
    span <- pmax(apply(xy, 2L, function(z) max(z) - min(z)), 1)
    xy <- sweep(sweep(xy, 2L, apply(xy, 2L, min), "-"), 2L, span, "/")
    distance <- as.matrix(stats::dist(xy))
    for (j in seq_len(nrow(pairs))) {
      a_name <- pairs[j, 1L]; b_name <- pairs[j, 2L]
      a <- as.logical(d[[a_name]]); b <- as.logical(d[[b_name]])
      if (sum(a) < 3L || sum(b) < 3L) next
      observed <- vapply(radii, function(r) .cross_pair_count(distance, a, b, r), numeric(1))
      null <- replicate(permutations, {
        pa <- sample(a); pb <- sample(b)
        vapply(radii, function(r) .cross_pair_count(distance, pa, pb, r), numeric(1))
      })
      if (is.null(dim(null))) null <- matrix(null, nrow = length(radii))
      expected <- rowMeans(null)
      for (r in seq_along(radii)) {
        k <- k + 1L
        p <- (1 + sum(abs(null[r, ] - expected[r]) >=
                        abs(observed[r] - expected[r]))) / (permutations + 1)
        rows[[k]] <- data.frame(patient = d$patient[1L], image = d$image[1L],
                                group = d$group[1L], pair = paste(a_name, b_name, sep = "--"),
                                radius_fraction = radii[r], observed_pairs = observed[r],
                                random_label_expected = expected[r],
                                log2_enrichment = log2((observed[r] + .5) / (expected[r] + .5)),
                                random_label_p = p, stringsAsFactors = FALSE)
      }
    }
  }
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
