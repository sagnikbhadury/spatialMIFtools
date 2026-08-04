#' Run the core patient-level spatial mIF workflow
#'
#' @param cells Physical-cell table.
#' @param phenotypes Binary phenotype columns.
#' @param group_levels Ordered c(reference, comparison) groups.
#' @param epithelial Epithelial phenotype used for zoning.
#' @param nx,ny Grid dimensions.
#' @param min_cells Minimum cells per retained bin.
#' @param sigma KDE smoothing width in bins.
#' @param ridge Ridge precision penalty.
#' @param minimum_network_bins Minimum bins per image network.
#' @param seed Random seed.
#' @return A named list of grids, zones, neighborhoods, networks, and contrasts.
#' @export
run_spatial_mif_pipeline <- function(cells, phenotypes, group_levels = NULL,
                                     epithelial = "Epithelial",
                                     nx = 12L, ny = 9L, min_cells = 5L,
                                     sigma = 1.25, ridge = 0.15,
                                     minimum_network_bins = 12L,
                                     seed = 20260803L) {
  validate_mif_cells(cells, phenotypes)
  grids <- bin_mif_cells(cells, phenotypes, nx, ny, min_cells)
  zones <- add_epithelial_zones(grids, epithelial, sigma, nx, ny)
  neighborhoods <- fit_spatial_neighborhoods(zones, phenotypes, seed = seed)
  assigned <- assign_neighborhoods(neighborhoods, zones)
  networks <- estimate_ridge_networks(assigned, phenotypes, ridge,
                                      minimum_network_bins, FALSE, FALSE)
  adjusted <- estimate_ridge_networks(assigned, phenotypes, ridge,
                                      minimum_network_bins, TRUE, FALSE)
  contrasts <- patient_edge_contrasts(networks, group_levels)
  adjusted_contrasts <- patient_edge_contrasts(adjusted, group_levels)
  structure(list(grids = grids, zones = zones, neighborhoods = neighborhoods,
                 assigned_grids = assigned, networks = networks,
                 adjusted_networks = adjusted, contrasts = contrasts,
                 adjusted_contrasts = adjusted_contrasts,
                 call = match.call()), class = "spatialmif_pipeline")
}

#' Optional adapter to SpatialMethodsWorkbench
#'
#' This function deliberately delegates to the installed workbench package and
#' does not vendor ISPAT or GP-GHS method code.
#'
#' @param data Input accepted by SpatialMethodsWorkbench::run_analysis().
#' @param method Workbench method identifier, for example "ISPAT" or "GPGHS".
#' @param ... Additional arguments forwarded to run_analysis().
#' @return The workbench result.
#' @export
spatial_methods_workbench_adapter <- function(data, method, ...) {
  if (!requireNamespace("SpatialMethodsWorkbench", quietly = TRUE)) {
    stop("Install SpatialMethodsWorkbench to use this adapter.", call. = FALSE)
  }
  do.call(SpatialMethodsWorkbench::run_analysis,
          c(list(data = data, method = method), list(...)))
}
