#' @title Posterior Distribution of Minimum X% Effective Dose
#' @param x output from `dreamer_mcmc()`.
#' @param ed a number between 0 and 100 indicating the ed% dose that is
#'   being sought.
#' @param probs a vector of quantiles to calculate on the posterior.
#' @param time the slice of time for which to calculate the posterior EDX dose.
#'   Applies to longitudinal models only.
#' @param lower the lower bound of the doses for calculating EDX.
#' @param upper the upper bound of the doses for calculating EDX.
#' @param greater if `TRUE`, higher values indicate better efficacy.  If
#'   `FALSE`, lower responses indicate better efficacy.
#' @param small_bound the minimum (`greater = TRUE`) or maximum
#'   (`greater = FALSE`) bound of the response.
#' @param index a vector indicating which MCMC samples to use in the
#'   calculation.  If `NULL` (default), all MCMC samples are used.
#' @param return_samples logical indicating if the posterior samples should be
#'   returned.
#' @param ... additional arguments for specific methods.
#' @details The minimum X% effective dose is the dose that has X% of the
#'   largest effect for doses between `lower` and `upper`.  When `greater`
#'   is `TRUE`, larger positive responses are considered more effective and
#'   vice versa.  The X% response is calculated as `small_bound` +
#'   `ed` / 100 * (max_effect - `small_bound`) where "max_effect" is the
#'   maximum response for doses between `lower` and `upper`.  The X% effective
#'   dose is the smallest dose which has X% response within the dose range.
#'   It is possible that for some MCMC samples, an X% effective dose may
#'   not exist, so probabilities are not guaranteed to sum to one.
#' @return Posterior quantities and samples (if applicable),
#'   generally in the form of a list.  The `pr_edx_exists` column gives the
#'   posterior probability that an EDX% effect exists.
#' @example man/examples/ex-post_medx.R
#' @export
post_medx <- function(
  x,
  ed,
  probs,
  time,
  lower,
  upper,
  greater,
  small_bound,
  return_samples,
  ...
) {
  UseMethod("post_medx", x)
}

#' @rdname post_medx
#' @export
post_medx.dreamer_bma <- function(
  x,
  ed,
  probs = c(.025, .975),
  time = NULL,
  lower = min(x$doses),
  upper = max(x$doses),
  greater = TRUE,
  small_bound = 0,
  return_samples = FALSE,
  ...
) {
  assert_no_dots("post_medx.dreamer_bma", ...)
  time <- get_time(x, time)
  mcmc_index <- vapply(
    x,
    function(y) any(grepl("mcmc_", class(y))),
    logical(1)
  ) %>% which()
  model_index <- attr(x, "model_index")
  samps <- purrr::map_df(
    unique(model_index),
    function(i) {
      index <- which(model_index == i)
      post_medx(
        x = x[[mcmc_index[i]]],
        ed = ed,
        probs = probs,
        time = time,
        lower = lower,
        upper = upper,
        greater = greater,
        small_bound = small_bound,
        index = index,
        return_samples = TRUE
      )$samps
    }
  )
  output <- summarize_medx(samps = samps, probs = probs)
  return_stats_samps(output, samps, return_samples)
}

#' @rdname post_medx
#' @export
post_medx.dreamer_mcmc <- function(
  x,
  ed,
  probs = c(.025, .975),
  time = NULL,
  lower = min(attr(x, "doses")),
  upper = max(attr(x, "doses")),
  greater = TRUE,
  small_bound = 0,
  return_samples = FALSE,
  index = 1:(nrow(x[[1]]) * length(x)),
  ...
) {
  assert_no_dots("post_medx.dreamer_mcmc", ...)
  time <- get_time(x, time)
  extremes <- get_extreme(
    x = x,
    time = time,
    greater = greater,
    lower = lower,
    upper = upper,
    index = index
  )
  samps <- purrr::map_df(
    ed,
    calc_post_medx_samps,
    effect100 = extremes$extreme_responses,
    extreme_dose = extremes$doses,
    x = x,
    index = index,
    lower = lower,
    upper = upper,
    greater = greater,
    time = time,
    small_bound = small_bound
  )
  output <- summarize_medx(samps = samps, probs = probs)
  return_stats_samps(output, samps, return_samples)
}

calc_post_medx_samps <- function(
  ed,
  effect100,
  extreme_dose,
  small_bound,
  x,
  index,
  lower,
  upper,
  greater,
  time
) {
  effect_x <- (ed / 100) * (effect100 - small_bound) + small_bound
  post_doses <- get_dose(
    x = x,
    time = time,
    response = effect_x,
    index = index,
    lower = lower,
    upper = upper
  )
  post_doses[post_doses > extreme_dose] <- NA
  return(
    dplyr::tibble(
      ed = ed,
      dose = post_doses
    )
  )
}

post_medx.dreamer_mcmc_independent <- function(...) { #nolint
  rlang::abort(
    "post_medx() not supported for independent models.",
    class = "dreamer"
  )
}

post_medx.dreamer_mcmc_independent_binary <- function(...) { #nolint
  rlang::abort(
    "post_medx() not supported for independent models.",
    class = "dreamer"
  )
}

summarize_medx <- function(samps, probs) {
  output <- dplyr::group_by(samps, .data$ed) %>%
    dplyr::summarize(
      pr_edx_exists = mean(!is.na(.data$dose)),
      mean = mean(.data$dose, na.rm = TRUE)
    )
  output2 <- purrr::map(
    probs,
    function(xx, samps) {
      dplyr::group_by(samps, .data$ed) %>%
        dplyr::summarize(
          !!(paste0(100 * xx, "%")) :=
            quantile(.data$dose, probs = xx, na.rm = TRUE)
        )
    },
    samps = samps
  ) %>%
    { #nolint
      Reduce(
        function(x, y, ...) merge(x, y, by = "ed"),
        .
      )
    }
  output <- merge(output, output2, by = "ed")
  return(output)
}

return_stats_samps <- function(stats, samps, return_samples) {
  if (return_samples) {
    return(list(stats = stats, samps = samps))
  } else {
    return(list(stats = stats))
  }
}
