#' @title Probability of minimum X\% effective dose
#' @description Calculate the probability a dose being the smallest dose
#'   that has at least X\% of the maximum efficacy.
#' @param x output from a call to `dreamer_mcmc()`.
#' @param doses the doses for which pr(minimum effective X\% dose) is to be
#'   calculated.
#' @param ed a number between 0 and 100 indicating the ed\% dose that is
#'   being sought.
#' @param greater if `TRUE`, higher responses indicate better efficacy.  If
#'   `FALSE`, lower responses indicate better efficacy.`
#' @param small_bound the lower (upper) bound of the response variable
#'   when `greater = TRUE` (`FALSE`).  This is used to calculate the
#'   `ed`\% effect as `ed / 100 * (effect_100 - small_bound) + small_bound`.
#' @param time the time (scalar) at which the Pr(MEDX) should be calculated.
#' @details Obtaining the probability of a particular does being the
#'   minimum efficacious dose achieving `ed`\% efficacy is dependent on
#'   the doses specified.
#'
#'   For a given MCMC sample of parameters, the 100\% efficacy value is defined
#'   as the highest efficacy of the doses specified.  For each posterior draw
#'   of MCMC parameters, the minimum `ed`\% efficacious dose is defined as the
#'   lowest dose what has at least `ed`\% efficacy relative to the 100\%
#'   efficacy value.
#'
#'   The `ed`\% effect is calculated as
#'   `ed / 100 * (effect_100 - small_bound) + small_bound` where `effect_100`
#'   is the largest mean response among `doses` for a given MCMC iteration.
#' @return A data frame with the following columns:
#'   - `dose`: numeric dose levels.
#'   - `prob`: Prob(EDX | data) for each dose. Note: these probabilities do
#'     not necessarily sum to 1 because the EDX may not exist. In fact,
#'     Pr(EDX does not exist | data) = `1 - sum(prob)`.
#' @example man/examples/ex-pr_medx.R
#' @export
pr_medx <- function(
  x,
  doses = attr(x, "doses"),
  ed,
  greater = TRUE,
  small_bound = 0,
  time = NULL
) {
  time <- get_time(x, time)
  if (greater) {
    ie <- get(">")
    mm <- max
  } else {
    ie <- get("<")
    mm <- min
  }
  samps <- posterior(
    x = x,
    doses = doses,
    return_samples = TRUE,
    return_stats = FALSE,
    times = time,
  )$samps %>%
    dplyr::mutate(
      dose_index = vapply(
        .data$dose,
        function(xx) which(xx == !!doses),
        integer(1)
      )
    ) %>%
    dplyr::group_by(.data$iter) %>%
    dplyr::mutate(
      effect_100 = mm(.data$mean_response)
    )
  check_small_bound(small_bound, nrow(samps))
  samps2 <- purrr::map_df(
    ed,
    function(ed) {
      samps %>%
        dplyr::mutate(
          effect_X = !!ed / 100 * (.data$effect_100 - !!small_bound) +
            !!small_bound,
          effective_dose = ie(.data$mean_response, .data$effect_X)
        ) %>%
        dplyr::group_by(.data$iter) %>%
        dplyr::summarize(
          med_index = new_min(.data$dose_index[.data$effective_dose])
        ) %>%
        dplyr::mutate(ed = !!ed) %>%
        dplyr::ungroup()
    }
  )
  out <- samps2 %>%
    dplyr::group_by(ed) %>%
    dplyr::summarize(data = list({
      out <- tabulate(.data$med_index, nbins = length(!!doses))
      names(out) <- !!doses
      out / length(.data$med_index)
    })
    ) %>%
    tidyr::unnest_wider(col = .data$data) %>%
    tidyr::pivot_longer(cols = -.data$ed, names_to = "dose", values_to = "prob")
  if (!is.null(time)) {
    out <- out %>%
      dplyr::mutate(time = !!time)
  }
  return(out)
}

check_small_bound <- function(small_bound, n_mcmc) {
  if ((length(small_bound) != 1) & (length(small_bound) != n_mcmc)) {
    stop("small_bound must have length 1 or ", n_mcmc, ".", call. = FALSE)
  }
}
