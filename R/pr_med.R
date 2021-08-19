#' @title Pr(minimum efficacious dose)
#' @description Calculates the posterior probability that each specified doses
#'   are the minimum effective dose in the set; i.e. the smallest
#'   dose that has a clinically significant difference (CSD).
#' @param x output from a call to `dreamer_mcmc()`.
#' @param doses the doses for which pr(MED) is to be calculated.
#' @param csd the treatment effect that is clinically relevant.
#' @param reference_dose a single dose that is used as the reference when
#'   defining the MED relative to a dose (rather than in absolute terms).  When
#'   `reference_dose` is specified, this function calculates the
#'   posterior probability
#'   that each dose is the smallest dose such that
#'   (effect_dose - effect_reference_dose > CSD).
#' @param greater if `TRUE`, higher responses indicate better efficacy.  If
#'   `FALSE`, lower responses indicate better efficacy.`
#' @param time the time (scalar) at which the Pr(MED) should be calculated.
#'   Applies only to longitudinal models.
#' @return A tibble listing each dose and the posterior probability that 
#'   each dose is the minimum efficacious dose.
#' @example man/examples/ex-pr_med.R
#' @export
pr_med <- function( #nolint
  x,
  doses = attr(x, "doses"),
  csd = NULL,
  reference_dose = NULL,
  greater = TRUE,
  time = NULL
) {
  assert_doses(doses)
  assert_reference_dose(reference_dose)
  time <- get_time(x, time)
  samps <- posterior(
    x = x,
    doses = doses,
    reference_dose = reference_dose,
    return_samples = TRUE,
    return_stats = FALSE,
    times = time,
  )$samps
  samps %>%
    get_med(csd, greater) %>%
    dplyr::mutate(n = n()) %>%
    dplyr::group_by(n) %>%
    dplyr::count(dose = .data$med, name = "prob") %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(.data$dose)) %>%
    dplyr::mutate(prob = .data$prob / .data$n) %>%
    tidyr::complete(dose = !!doses, fill = list(prob = 0)) %>%
    dplyr::select(.data$dose, .data$prob) %>%
    add_cols(reference_dose, time) %>%
    dplyr::select(.data$dose, everything(), -.data$prob, .data$prob)
}

add_cols <- function(x, reference_dose, time) {
  if (!is.null(reference_dose)) {
    x <- x %>%
      dplyr::mutate(reference_dose = !!reference_dose)
  }
  if (!is.null(time)) {
    x <- x %>%
      dplyr::mutate(time = !!time)
  }
  x
}

add_efficacious <- function(x, threshold, greater) {
  ie <- get_inequality_func(greater)
  dplyr::mutate(
    x,
    efficacious = ie(.data$mean_response, !!threshold)
  )
}

get_med <- function(x, threshold, greater) {
  add_efficacious(x, threshold, greater) %>%
    dplyr::group_by(.data$iter) %>%
    dplyr::summarize(
      med = new_min(.data$dose[.data$efficacious]),
      .groups = "drop"
    )
}

get_inequality_func <- function(greater) {
  if (greater) {
    ie <- get(">")
  } else {
    ie <- get("<")
  }
  ie
}

get_time <- function(x, time, max_length = 1) {
  if (is.null(time) & !is.null(attr(x, "times"))) {
    time <- max(attr(x, "times"))
  } else if (is.null(time) & !is.null(x$times)) {
    time <- max(x$times)
  } else if (!is.null(time)) {
    if (length(time) > max_length) {
      rlang::abort(
        paste0("argument 'time' must have length <= ", max_length, "."),
        class = "dreamer"
      )
    }
  }
  return(time)
}

new_min <- function(x) {
  if (length(x) > 0) {
    return(min(x))
  } else{
    return(NA_real_)
  }
}
