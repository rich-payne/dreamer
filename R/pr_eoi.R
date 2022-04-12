#' @title Calculate Probability of Meeting Effect of Interest (EOI)
#' @description Calculate Pr(effect_dose - effect_reference_dose > EOI | data)
#'  or Pr(effect_dose > EOI | data).
#' @param x output from a call to `dreamer_mcmc()`.
#' @param eoi a vector of the effects of interest (EOI) in the probability
#'   function.
#' @param dose a vector of the doses for which to calculate the
#'   posterior probabilities.
#' @param reference_dose a vector of doses for relative effects of interest.
#' @param time the time at which to calculate the posterior quantity.  Defaults
#'   to the latest timepoint.  Applies to longitudinal models only.
#' @return A tibble listing the doses, times, and
#'   Pr(effect_dose - effect_reference_dose > eoi) if `reference_dose`
#'   is specified; otherwise, Pr(effect_dose > eoi).
#' @example man/examples/ex-pr_eoi.R
#' @export
pr_eoi <- function(x, eoi, dose, reference_dose = NULL, time = NULL) {
  time <- get_time(x, time)
  grid <- get_eoi_grid(eoi, dose, reference_dose)
  output <- purrr::pmap_dfr(
    grid,
    pr_eoi_impl,
    time = time,
    x = x
  )
  output %>%
    dplyr::select(
      any_of(c("eoi", "time")), any_of(c("dose", "reference_dose", "prob"))
    )
}

pr_eoi_impl <- function(x, eoi, dose, reference_dose = NULL, time) {
  samps <- posterior(
    x = x,
    doses = dose,
    reference_dose = reference_dose,
    times = time,
    return_samples = TRUE,
    return_stats = FALSE
  )$samps
  tibble::tibble(
    eoi = eoi,
    dose = dose,
    time = !!time,
    reference_dose = !!reference_dose,
    prob = mean(samps$mean_response > eoi)
  )
}

get_eoi_grid <- function(eoi, dose, reference_dose) {
  n_eoi <- length(eoi)
  n_dose <- length(dose)
  n_ref <- length(reference_dose)
  lens <- c(n_eoi, n_dose, n_ref)
  ind <- lens > 1
  lens_gt_1 <- lens[ind]
  n_gt_1 <- sum(ind)
  if (n_gt_1 == 0) {
    same_lengths <- TRUE
  } else if (n_gt_1 > 0) {
    same_lengths <- all(lens_gt_1 == lens_gt_1[1])
  }
  if (n_gt_1 == 0 || same_lengths) {
    grid <- tibble::tibble(eoi, dose, reference_dose = !!reference_dose)
  } else if (n_gt_1 > 1 && !same_lengths) {
    grid <- tidyr::expand_grid(eoi, dose, reference_dose = !!reference_dose)
  }
  grid
}
