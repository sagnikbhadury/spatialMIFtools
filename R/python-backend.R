#' Test availability of the optional Python backend
#' @param python Optional Python executable passed to reticulate.
#' @return TRUE when reticulate and required modules are available.
#' @export
python_backend_available <- function(python = NULL) {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(FALSE)
  if (!is.null(python)) reticulate::use_python(python, required = FALSE)
  all(vapply(c("numpy", "pandas", "scipy"), reticulate::py_module_available,
             logical(1)))
}

#' Load the bundled Python analysis backend
#'
#' The backend contains computational functions extracted from the manuscript's
#' validated Python pipeline. R data frames and matrices are converted by
#' reticulate.
#'
#' @param python Optional Python executable.
#' @param required Fail when Python dependencies are unavailable.
#' @return A reticulate Python module.
#' @export
spatial_mif_python <- function(python = NULL, required = TRUE) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    if (required) stop("Install the suggested package 'reticulate'.", call. = FALSE)
    return(NULL)
  }
  if (!is.null(python)) reticulate::use_python(python, required = required)
  path <- system.file("python", package = "spatialMIFtools")
  if (!nzchar(path)) {
    candidate <- file.path(getwd(), "inst", "python")
    if (dir.exists(candidate)) path <- normalizePath(candidate)
  }
  if (!nzchar(path) || !dir.exists(path)) stop("Bundled Python backend not found.", call. = FALSE)
  reticulate::import_from_path("spatial_mif_backend", path = path,
                               convert = TRUE, delay_load = !required)
}
