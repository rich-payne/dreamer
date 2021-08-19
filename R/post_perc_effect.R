#' @title Calculate Posterior of a Dose's Percentage Effect
#' @description Given a dose, the "percentage effect" is defined as
#'   (effect of the given dose - small_bound) / (maximum effect in dose range - 
#'   small_bound).  This function returns the posterior statistics and/or
#'   samples of this effect.
#' @param x output from a call to `dreamer_mcmc()`, or
#'   the MCMC samples from a single model of output from
#'   a `dreamer_mcmc()` call.
#' @param dose the dose at which to calculate the posterior percentage
#'   effect.
#' @param lower the lower bound of the dose range under consideration.
#' @param upper the upper bound of the dose range under consideration.
#' @param greater logical indicating if the response is desired to
#'   be increasing (`TRUE`) or decreasing (`FALSE`).
#' @param small_bound the lower (if `greater = TRUE`) or upper
#'   (if `greater = FALSE`) bound that the effect is expected to
#'   take.
#' @param index an index on which MCMC samples should be used.  Generally
#'   the user should not specify anything for this argument as `dreamer`
#'   will handle this automatically.
#' @param probs a vector of quantiles to calculate on the posterior.
#' @param time the slice of time for which to calculate the posterior
#'   percentage effect.  Applies to longitudinal models only.
#' @param return_samples logical indicating if the posterior samples should
#'   be returned.
#' @return A named list with the following components:
#'   * stats: a tibble listing the dose, time (where relevant),
#'     probability a percentage effect exists, the average percentage effect,
#'     and the specified quantiles of the percentage effect.
#'   * samps: a tibble with the posterior samples for each dose/time
#'     combination.
#' @example man/examples/ex-post_perc_effect.R
#' @export
post_perc_effect <- function(
  x,
  dose,
  probs,
  time,
  lower,
  upper,
  greater,
  small_bound,
  index,
  return_samples
) {
  UseMethod("post_perc_effect", x)
}

#' @rdname post_perc_effect
#' @export
post_perc_effect.dreamer_bma <- function(
  x,
  dose,
  probs = c(.025, .975),
  time = NULL,
  lower = min(x$doses),
  upper = max(x$doses),
  greater = TRUE,
  small_bound = 0,
  index = NA,
  return_samples = FALSE
) {
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
      post_perc_effect(
        x = x[[mcmc_index[i]]],
        dose = dose,
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
  output <- summarize_post_perc(samps = samps, probs = probs)
  return_stats_samps(output, samps, return_samples)
}

#' @rdname post_perc_effect
#' @export
post_perc_effect.dreamer <- function(
  x,
  dose,
  probs = c(.025, .975),
  time = NULL,
  lower = min(attr(x, "doses")),
  upper = max(attr(x, "doses")),
  greater = TRUE,
  small_bound = 0,
  index = 1:(nrow(x[[1]]) * length(x)),
  return_samples = FALSE
) {
  time <- get_time(x, time)
  check_bounds(dose, lower, upper)
  check_small_bound2(small_bound, index)
  extreme_responses <- get_extreme(
    x = x,
    time = time,
    greater = greater,
    lower = lower,
    upper = upper,
    index = index
  )$extreme_responses
  samps <- purrr::map_df(
    dose,
    function(
      dose,
      x,
      index,
      small_bound,
      extreme_responses,
      greater,
      time
    ) {
      dose_response <- dreamer_mean(
        x = x,
        dose = dose,
        time = time,
        index = index
      )
      post_perc <- (dose_response - small_bound) /
        (extreme_responses - small_bound)
      post_perc[post_perc < 0 | post_perc > 1] <- NA
      return(dplyr::tibble(dose = dose, perc_effect = post_perc))
    },
    x = x,
    index = index,
    small_bound = small_bound,
    extreme_responses = extreme_responses,
    greater = greater,
    time = time
  )
  output <- summarize_post_perc(samps = samps, probs = probs)
  return_stats_samps(output, samps, return_samples)
}

summarize_post_perc <- function(samps, probs) {
  output <- dplyr::group_by(samps, .data$dose) %>%
    dplyr::summarize(
      pr_perc_exists = mean(!is.na(.data$perc_effect)),
      mean = mean(.data$perc_effect, na.rm = TRUE)
    )
  output2 <- purrr::map(
    probs,
    function(xx, samps) {
      nme <- sprintf("%.2f%%", 100 * xx)
      dplyr::group_by(samps, .data$dose) %>%
        dplyr::summarize(
          !!nme := quantile(
            .data$perc_effect,
            probs = xx,
            na.rm = TRUE,
            names = FALSE
          )
        )
    },
    samps = samps
  ) %>%
    { #nolint
      Reduce(
        function(x, y, ...) merge(x, y, by = "dose"),
        .
      )
    }
  output <- merge(output, output2, by = "dose") %>%
    dplyr::as_tibble()
  return(output)
}

check_bounds <- function(dose, lower, upper) {
  if (any(dose < lower) || any(dose > upper)) {
    stop("dose must be between lower and upper")
  }
}

check_small_bound2 <- function(small_bound, index) {
  if (!((length(small_bound)) == 1 || length(small_bound) == length(index))) {
    stop("Small_bound must have length 1 or length(index).", call. = FALSE)
  }
}
