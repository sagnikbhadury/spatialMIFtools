.required_cell_columns <- c("patient", "image", "group", "x", "y")

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

.check_columns <- function(x, columns, object = deparse(substitute(x))) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) {
    stop(object, " is missing required columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  invisible(TRUE)
}

.assert_positive_integer <- function(x, name) {
  if (length(x) != 1L || !is.finite(x) || x < 1 || x != as.integer(x)) {
    stop(name, " must be one positive integer.", call. = FALSE)
  }
  as.integer(x)
}

.safe_scale <- function(x) {
  x <- as.matrix(x)
  center <- colMeans(x, na.rm = TRUE)
  scale <- apply(x, 2L, stats::sd, na.rm = TRUE)
  active <- is.finite(scale) & scale > 0
  z <- matrix(0, nrow(x), ncol(x), dimnames = dimnames(x))
  if (any(active)) {
    z[, active] <- sweep(sweep(x[, active, drop = FALSE], 2L,
                               center[active], "-"), 2L, scale[active], "/")
  }
  list(x = z, center = center, scale = scale, active = active)
}

.apply_scaler <- function(x, scaler) {
  x <- as.matrix(x)
  out <- matrix(0, nrow(x), ncol(x), dimnames = dimnames(x))
  active <- scaler$active & colnames(x) %in% names(scaler$center)
  if (any(active)) {
    cols <- colnames(x)[active]
    out[, active] <- sweep(sweep(x[, active, drop = FALSE], 2L,
                                 scaler$center[cols], "-"), 2L,
                          scaler$scale[cols], "/")
  }
  out
}

.auc <- function(y, score) {
  y <- as.integer(y)
  ok <- is.finite(score) & y %in% c(0L, 1L)
  y <- y[ok]
  score <- score[ok]
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  if (!n1 || !n0) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

.ridge_logistic_fit <- function(x, y, lambda = 4) {
  x <- as.matrix(x)
  y <- as.numeric(y)
  design <- cbind("(Intercept)" = 1, x)
  objective <- function(beta) {
    eta <- drop(design %*% beta)
    nll <- sum(log1p(exp(-abs(eta))) + pmax(eta, 0) - y * eta)
    nll + 0.5 * lambda * sum(beta[-1L]^2)
  }
  gradient <- function(beta) {
    eta <- drop(design %*% beta)
    p <- stats::plogis(eta)
    grad <- drop(crossprod(design, p - y))
    grad[-1L] <- grad[-1L] + lambda * beta[-1L]
    grad
  }
  fit <- stats::optim(rep(0, ncol(design)), objective, gradient,
                      method = "BFGS", control = list(maxit = 1000))
  structure(list(coefficients = fit$par, convergence = fit$convergence,
                 columns = colnames(x)), class = "spatialmif_ridge_logistic")
}

.ridge_logistic_predict <- function(object, newx) {
  newx <- as.matrix(newx)[, object$columns, drop = FALSE]
  stats::plogis(drop(cbind(1, newx) %*% object$coefficients))
}

.pair_names <- function(phenotypes) {
  pairs <- utils::combn(phenotypes, 2L)
  data.frame(source = pairs[1L, ], target = pairs[2L, ],
             edge = paste(pairs[1L, ], pairs[2L, ], sep = "--"),
             stringsAsFactors = FALSE)
}

.image_split <- function(x) split(x, interaction(x$patient, x$image, drop = TRUE))
