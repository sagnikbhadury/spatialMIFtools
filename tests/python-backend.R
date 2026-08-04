library(spatialMIFtools)

python <- Sys.getenv("SPATIALMIFTOOLS_PYTHON", unset = "")
if (nzchar(python) && requireNamespace("reticulate", quietly = TRUE)) {
  py <- spatial_mif_python(python = python)
  set.seed(44)
  x <- matrix(rnorm(500), ncol = 5)
  r_value <- ridge_partial_correlation(x, ridge = .15)
  py_value <- py$partial_correlation(x, ridge = .15)
  stopifnot(max(abs(r_value - py_value), na.rm = TRUE) < 1e-10)

  phenotypes <- c("Epithelial", "APC", "Treg", "CTL", "HelperT")
  cells <- simulate_mif_cells(4, 1, 100, seed = 45)
  py_grid <- py$grid_image(cells, phenotypes, nx = 8L, ny = 6L,
                           minimum_cells = 2L)
  stopifnot(nrow(py_grid) > 0,
            all(paste0("density__", phenotypes) %in% names(py_grid)))
}
