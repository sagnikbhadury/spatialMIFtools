#' Patient-clustered ordered zone interaction
#'
#' @param data Long-form image-zone summary data.
#' @param value Numeric outcome column.
#' @param patient,image,group,zone Column names.
#' @param group_levels Ordered c(reference, comparison) groups.
#' @param corstr GEE working correlation structure.
#' @return The group-by-zone coefficient and robust uncertainty.
#' @export
ordered_zone_gee <- function(data, value, patient = "patient", image = "image",
                             group = "group", zone = "tumor_zone",
                             group_levels = NULL, corstr = "exchangeable") {
  if (!requireNamespace("geepack", quietly = TRUE)) {
    stop("ordered_zone_gee() requires the suggested package 'geepack'.",
         call. = FALSE)
  }
  .check_columns(data, c(value, patient, image, group, zone), "data")
  d <- data[stats::complete.cases(data[c(value, patient, group, zone)]), ]
  zone_map <- c(Q1_low = 1, Q2 = 2, Q3 = 3, Q4_high = 4)
  if (is.numeric(d[[zone]])) {
    d$.zone_score <- as.numeric(d[[zone]])
  } else {
    d$.zone_score <- unname(zone_map[as.character(d[[zone]])])
  }
  d <- d[is.finite(d$.zone_score), ]
  groups <- unique(as.character(d[[group]]))
  group_levels <- group_levels %||% groups[seq_len(min(2L, length(groups)))]
  if (length(group_levels) != 2L) stop("Exactly two groups are required.", call. = FALSE)
  d$.group_indicator <- as.integer(as.character(d[[group]]) == group_levels[2L])
  d$.patient_id <- as.factor(d[[patient]])
  d$.outcome <- as.numeric(d[[value]])
  fit <- geepack::geeglm(.outcome ~ .group_indicator * .zone_score,
                         id = d$.patient_id, data = d, family = stats::gaussian,
                         corstr = corstr, std.err = "san.se")
  tab <- summary(fit)$coefficients
  term <- ".group_indicator:.zone_score"
  if (!term %in% rownames(tab)) stop("Interaction coefficient was not estimable.", call. = FALSE)
  estimate <- unname(tab[term, "Estimate"])
  se_name <- intersect(c("Std.err", "Std. Error"), colnames(tab))[1L]
  p_name <- grep("Pr\\(", colnames(tab), value = TRUE)[1L]
  se <- unname(tab[term, se_name])
  data.frame(reference = group_levels[1L], comparison = group_levels[2L],
             estimate = estimate, std_error = se,
             conf_low = estimate - 1.96 * se, conf_high = estimate + 1.96 * se,
             p_value = unname(tab[term, p_name]), n_patients = length(unique(d$.patient_id)),
             stringsAsFactors = FALSE)
}

#' Patient-label permutation for a vector of statistics
#'
#' @param data Analysis data containing patient and group columns.
#' @param statistic Function taking a data frame and returning a named numeric vector.
#' @param patient,group Column names.
#' @param permutations Number of patient-label permutations.
#' @param target Optional name of the prespecified target statistic.
#' @param seed Random seed.
#' @return Observed statistics, permutation matrix, and target/familywise empirical p-values.
#' @export
patient_label_permutation <- function(data, statistic, patient = "patient",
                                      group = "group", permutations = 999L,
                                      target = NULL, seed = 20260803L) {
  .check_columns(data, c(patient, group), "data")
  if (!is.function(statistic)) stop("statistic must be a function.", call. = FALSE)
  mapping <- unique(data[c(patient, group)])
  if (anyDuplicated(mapping[[patient]])) {
    stop("Each patient must map to one group.", call. = FALSE)
  }
  observed <- statistic(data)
  if (!is.numeric(observed) || is.null(names(observed))) {
    stop("statistic(data) must return a named numeric vector.", call. = FALSE)
  }
  permutations <- .assert_positive_integer(permutations, "permutations")
  set.seed(seed)
  null <- matrix(NA_real_, permutations, length(observed),
                 dimnames = list(NULL, names(observed)))
  for (b in seq_len(permutations)) {
    perm_map <- mapping
    perm_map[[group]] <- sample(mapping[[group]])
    d <- data
    d[[group]] <- perm_map[[group]][match(d[[patient]], perm_map[[patient]])]
    value <- statistic(d)
    null[b, names(value)] <- value
  }
  result <- list(observed = observed, permutation = null)
  if (!is.null(target)) {
    if (!target %in% names(observed)) stop("target is not returned by statistic.", call. = FALSE)
    result$target_p <- (1 + sum(abs(null[, target]) >= abs(observed[target]), na.rm = TRUE)) /
      (permutations + 1)
    maxima <- apply(abs(null), 1L, max, na.rm = TRUE)
    result$familywise_p <- (1 + sum(maxima >= abs(observed[target]))) /
      (permutations + 1)
  }
  class(result) <- "spatialmif_permutation"
  result
}
