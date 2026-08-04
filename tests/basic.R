library(spatialMIFtools)

phenotypes <- c("Epithelial", "APC", "Treg", "CTL", "HelperT")
cells <- simulate_mif_cells(8, 2, 120, seed = 101)
stopifnot(nrow(cells) > 0)

grids <- bin_mif_cells(cells, phenotypes, nx = 8, ny = 6, min_cells = 2)
stopifnot(nrow(grids) > 0, all(paste0("density__", phenotypes) %in% names(grids)))

zones <- add_epithelial_zones(grids, sigma = 1)
stopifnot(all(c("epithelial_kde", "tumor_zone") %in% names(zones)))

pc <- ridge_partial_correlation(as.matrix(zones[paste0("density__", phenotypes)]))
stopifnot(all(dim(pc) == c(5, 5)), all(diag(pc) == 1))

networks <- estimate_ridge_networks(zones, phenotypes, minimum_bins = 5)
stopifnot(nrow(networks) > 0)
contrasts <- patient_edge_contrasts(networks, c("Group1", "Group2"))
stopifnot(nrow(contrasts) == choose(length(phenotypes), 2))

nh <- fit_spatial_neighborhoods(zones, phenotypes, k_range = 3:4,
                                max_silhouette_bins = 500, seed = 102)
assigned <- assign_neighborhoods(nh, zones)
stopifnot("neighborhood" %in% names(assigned))

scores <- nested_lopo_prediction(zones, phenotypes,
                                 group_levels = c("Group1", "Group2"),
                                 feature_set = "neighborhood", k_range = 3:4,
                                 seed = 103)
stopifnot(nrow(scores) == 8, is.finite(attr(scores, "auc")))

zone_data <- expand.grid(patient = sprintf("P%02d", 1:12), image = 1:2,
                         tumor_zone = c("Q1_low", "Q2", "Q3", "Q4_high"),
                         stringsAsFactors = FALSE)
zone_data$group <- ifelse(as.integer(sub("P", "", zone_data$patient)) <= 6,
                          "Group1", "Group2")
zone_data$value <- as.numeric(factor(zone_data$tumor_zone,
                                     levels = c("Q1_low", "Q2", "Q3", "Q4_high"))) *
  ifelse(zone_data$group == "Group2", -.05, 0) + rnorm(nrow(zone_data), sd = .1)
if (requireNamespace("geepack", quietly = TRUE)) {
  gee <- ordered_zone_gee(zone_data, "value", group_levels = c("Group1", "Group2"))
  stopifnot(nrow(gee) == 1, is.finite(gee$estimate))
}
