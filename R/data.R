#' Validate a multiplex imaging cell table
#'
#' @param cells Data frame with patient, image, group, x, y, and phenotype columns.
#' @param phenotypes Character vector naming binary phenotype columns.
#' @param allow_overlap Whether one cell may be positive for multiple phenotypes.
#' @return The input invisibly, after validation.
#' @export
validate_mif_cells <- function(cells, phenotypes, allow_overlap = TRUE) {
  if (!is.data.frame(cells)) stop("cells must be a data frame.", call. = FALSE)
  if (!length(phenotypes) || anyDuplicated(phenotypes)) {
    stop("phenotypes must contain unique column names.", call. = FALSE)
  }
  .check_columns(cells, c(.required_cell_columns, phenotypes), "cells")
  if (any(!is.finite(cells$x)) || any(!is.finite(cells$y))) {
    stop("x and y must be finite for every row.", call. = FALSE)
  }
  if (anyNA(cells$patient) || anyNA(cells$image) || anyNA(cells$group)) {
    stop("patient, image, and group may not be missing.", call. = FALSE)
  }
  bad <- vapply(cells[phenotypes], function(z) {
    z <- as.numeric(z)
    anyNA(z) || any(!z %in% c(0, 1))
  }, logical(1))
  if (any(bad)) {
    stop("Phenotype columns must contain only non-missing 0/1 values: ",
         paste(names(bad)[bad], collapse = ", "), call. = FALSE)
  }
  if (!allow_overlap && any(rowSums(cells[phenotypes]) > 1)) {
    stop("Overlapping phenotype calls found while allow_overlap = FALSE.",
         call. = FALSE)
  }
  invisible(cells)
}

#' Simulate a small multiplex imaging cohort
#'
#' Generates physical-cell coordinates and five binary annotations for examples
#' and software tests. It is not intended to reproduce either manuscript cohort.
#'
#' @param n_patients Number of patients.
#' @param images_per_patient Number of images per patient.
#' @param cells_per_image Approximate cells per image.
#' @param seed Random seed.
#' @return A cell-level data frame.
#' @export
simulate_mif_cells <- function(n_patients = 12L, images_per_patient = 2L,
                               cells_per_image = 180L, seed = 20260803L) {
  n_patients <- .assert_positive_integer(n_patients, "n_patients")
  images_per_patient <- .assert_positive_integer(images_per_patient,
                                                  "images_per_patient")
  cells_per_image <- .assert_positive_integer(cells_per_image, "cells_per_image")
  set.seed(seed)
  phenotypes <- c("Epithelial", "APC", "Treg", "CTL", "HelperT")
  rows <- vector("list", n_patients * images_per_patient)
  k <- 0L
  for (p in seq_len(n_patients)) {
    group <- if (p <= floor(n_patients / 2)) "Group1" else "Group2"
    for (i in seq_len(images_per_patient)) {
      k <- k + 1L
      n <- stats::rpois(1L, cells_per_image)
      centers <- matrix(c(0.28, 0.30, 0.72, 0.66, 0.52, 0.48), ncol = 2,
                        byrow = TRUE)
      component <- sample(1:3, n, replace = TRUE,
                          prob = if (group == "Group1") c(.45, .30, .25) else c(.25, .48, .27))
      xy <- centers[component, , drop = FALSE] + matrix(stats::rnorm(2 * n, sd = .12), ncol = 2)
      xy[xy < 0] <- 0
      xy[xy > 1] <- 1
      xy <- xy * 1000
      epithelial_p <- stats::plogis((xy[, 1] - 480) / 110)
      helper_p <- stats::plogis((600 - xy[, 1]) / 180) * .45
      ctl_p <- stats::plogis((xy[, 2] - 430) / 150) * if (group == "Group2") .55 else .40
      apc_p <- .12 + .35 * (component == 1)
      treg_p <- rep(.08, n)
      calls <- cbind(
        Epithelial = stats::rbinom(n, 1, epithelial_p),
        APC = stats::rbinom(n, 1, pmin(apc_p, .9)),
        Treg = stats::rbinom(n, 1, treg_p),
        CTL = stats::rbinom(n, 1, pmin(ctl_p, .9)),
        HelperT = stats::rbinom(n, 1, pmin(helper_p, .9))
      )
      rows[[k]] <- data.frame(
        patient = sprintf("P%03d", p), image = sprintf("P%03d_I%02d", p, i),
        group = group, x = xy[, 1], y = xy[, 2], calls,
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  validate_mif_cells(out, phenotypes)
  out
}
