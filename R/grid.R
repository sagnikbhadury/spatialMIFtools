#' Aggregate physical cells to area-normalized spatial grids
#'
#' @param cells Validated cell-level data frame.
#' @param phenotypes Binary phenotype columns.
#' @param nx,ny Number of horizontal and vertical bins.
#' @param min_cells Minimum total cells required to retain a bin.
#' @param scale_factor Density scaling constant before log1p transformation.
#' @return One row per retained image bin.
#' @export
bin_mif_cells <- function(cells, phenotypes, nx = 12L, ny = 9L,
                          min_cells = 5L, scale_factor = 1e4) {
  validate_mif_cells(cells, phenotypes)
  nx <- .assert_positive_integer(nx, "nx")
  ny <- .assert_positive_integer(ny, "ny")
  min_cells <- .assert_positive_integer(min_cells, "min_cells")
  if (!is.finite(scale_factor) || scale_factor <= 0) {
    stop("scale_factor must be positive.", call. = FALSE)
  }
  pieces <- .image_split(cells)
  out <- lapply(pieces, function(d) {
    xr <- max(max(d$x) - min(d$x), 1)
    yr <- max(max(d$y) - min(d$y), 1)
    bx <- pmin(floor(nx * (d$x - min(d$x)) / xr), nx - 1L)
    by <- pmin(floor(ny * (d$y - min(d$y)) / yr), ny - 1L)
    key <- interaction(bx, by, drop = TRUE)
    ids <- split(seq_len(nrow(d)), key)
    area <- max((xr / nx) * (yr / ny), 1)
    rows <- lapply(ids, function(ix) {
      if (length(ix) < min_cells) return(NULL)
      counts <- colSums(d[ix, phenotypes, drop = FALSE])
      ans <- data.frame(
        patient = as.character(d$patient[ix[1L]]),
        image = as.character(d$image[ix[1L]]),
        group = as.character(d$group[ix[1L]]),
        bin_x = bx[ix[1L]], bin_y = by[ix[1L]],
        x = (bx[ix[1L]] + .5) / nx, y = (by[ix[1L]] + .5) / ny,
        n_cells = length(ix), image_cells = nrow(d), bin_area = area,
        stringsAsFactors = FALSE
      )
      if ("cohort" %in% names(d)) ans$cohort <- as.character(d$cohort[ix[1L]])
      for (f in phenotypes) {
        ans[[f]] <- unname(counts[f])
        ans[[paste0("density__", f)]] <- log1p(counts[f] / area * scale_factor)
        ans[[paste0("prop__", f)]] <- counts[f] / length(ix)
      }
      ans
    })
    do.call(rbind, rows)
  })
  out <- Filter(Negate(is.null), out)
  if (!length(out)) return(data.frame())
  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  attr(ans, "phenotypes") <- phenotypes
  attr(ans, "grid") <- c(nx = nx, ny = ny)
  ans
}

.gaussian_kernel <- function(sigma) {
  radius <- max(1L, ceiling(4 * sigma))
  x <- -radius:radius
  k <- exp(-0.5 * (x / sigma)^2)
  k / sum(k)
}

.smooth_matrix <- function(x, sigma) {
  k <- .gaussian_kernel(sigma)
  radius <- (length(k) - 1L) / 2L
  pad_apply <- function(mat, margin) {
    apply(mat, margin, function(v) {
      padded <- c(rep(v[1L], radius), v, rep(v[length(v)], radius))
      stats::filter(padded, k, sides = 2L)[seq_along(v) + radius]
    })
  }
  first <- t(pad_apply(x, 1L))
  pad_apply(first, 2L)
}

#' Add within-image epithelial-density KDE quartiles
#'
#' @param grids Output from [bin_mif_cells()].
#' @param epithelial Phenotype count column used to construct the surface.
#' @param sigma Gaussian standard deviation in grid-bin units.
#' @param nx,ny Grid dimensions. Defaults to attributes saved by [bin_mif_cells()].
#' @return The grid with epithelial_kde and tumor_zone columns.
#' @export
add_epithelial_zones <- function(grids, epithelial = "Epithelial", sigma = 1.25,
                                 nx = NULL, ny = NULL) {
  .check_columns(grids, c("patient", "image", "bin_x", "bin_y", epithelial), "grids")
  dims <- attr(grids, "grid")
  nx <- nx %||% unname(dims["nx"])
  ny <- ny %||% unname(dims["ny"])
  if (!length(nx) || !length(ny) || any(!is.finite(c(nx, ny)))) {
    nx <- max(grids$bin_x) + 1L
    ny <- max(grids$bin_y) + 1L
  }
  if (!is.finite(sigma) || sigma <= 0) stop("sigma must be positive.", call. = FALSE)
  pieces <- .image_split(grids)
  out <- lapply(pieces, function(d) {
    surface <- matrix(0, nrow = ny, ncol = nx)
    surface[cbind(d$bin_y + 1L, d$bin_x + 1L)] <- d[[epithelial]]
    smooth <- .smooth_matrix(surface, sigma)
    values <- smooth[cbind(d$bin_y + 1L, d$bin_x + 1L)]
    d$epithelial_kde <- values
    if (!length(unique(values)) || max(values) <= 0 || stats::sd(values) == 0) {
      d$tumor_zone <- "Q0_no_epithelial"
    } else {
      ranks <- rank(values, ties.method = "first")
      cutpoints <- stats::quantile(ranks, probs = 0:4 / 4, type = 1)
      cutpoints <- unique(cutpoints)
      if (length(cutpoints) < 5L) {
        idx <- pmin(4L, ceiling(4 * ranks / length(ranks)))
      } else {
        idx <- cut(ranks, breaks = cutpoints, include.lowest = TRUE,
                   labels = FALSE)
      }
      d$tumor_zone <- c("Q1_low", "Q2", "Q3", "Q4_high")[idx]
    }
    d
  })
  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  attr(ans, "phenotypes") <- attr(grids, "phenotypes")
  attr(ans, "grid") <- c(nx = nx, ny = ny)
  ans
}
