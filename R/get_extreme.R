# get the extremes (min/max) of a dose response over a dose range

get_extreme <- function(x, time, greater, lower, upper, index) {
  UseMethod("get_extreme", x)
}

#' @keywords internal
#' @export
get_extreme.default <- function(x, time, greater, lower, upper, index) {
  rlang::abort(
    paste0(
      "Class(es) not supported for get_extreme: ",
      paste0(class(x), collapse = ", "),
      "."
    ),
    class = "dreamer"
  )
}

get_extreme_monotonic <- function(
    x,
    time,
    greater,
    lower,
    upper,
    index
) {
  check_time(time)
  doses <- c(lower, upper)
  responses <- posterior(
    x,
    doses = doses,
    times = time,
    return_samples = TRUE,
    iter = index,
    return_stats = FALSE
  )$samps
  er <- pextreme(
    responses = responses,
    greater = greater
  )
  return(er)
}

#' @keywords internal
#' @export
get_extreme.dreamer_mcmc_linear <- get_extreme_monotonic # nolint

#' @keywords internal
#' @export
get_extreme.dreamer_mcmc_loglinear <- get_extreme_monotonic # nolint

#' @keywords internal
#' @export
get_extreme.dreamer_mcmc_emax <- get_extreme_monotonic

#' @keywords internal
#' @export
get_extreme.dreamer_mcmc_exp <- get_extreme_monotonic

#' @keywords internal
#' @export
get_extreme.dreamer_mcmc_linear_binary <- get_extreme_monotonic # nolint

#' @keywords internal
#' @export
get_extreme.dreamer_mcmc_loglinear_binary <- get_extreme_monotonic # nolint

#' @keywords internal
#' @export
get_extreme.dreamer_mcmc_emax_binary <- get_extreme_monotonic # nolint

#' @keywords internal
#' @export
get_extreme.dreamer_mcmc_exp_binary <- get_extreme_monotonic # nolint

#' @keywords internal
#' @export
get_extreme.dreamer_mcmc_independent <- function( #nolint
  x,
  time,
  greater,
  lower = min(attr(x, "doses")),
  upper = max(attr(x, "doses")),
  index
) {
  check_time(time)
  doses <- attr(x, "doses")
  doses <- doses[doses >= lower & doses <= upper]
  responses <- posterior(
    x,
    doses = doses,
    times = time,
    return_samples = TRUE,
    iter = index,
    return_stats = FALSE
  )$samps
  er <- pextreme(
    responses = responses,
    greater = greater
  )
  return(er)
}

#' @keywords internal
#' @export
get_extreme.dreamer_mcmc_independent_binary <- # nolint
  get_extreme.dreamer_mcmc_independent

#' @keywords internal
#' @export
get_extreme.dreamer_mcmc_quad <- function(
  x,
  time,
  greater,
  lower,
  upper,
  index
) {
  check_time(time)
  doses1 <- c(lower, upper)
  y <- subset_mcmc(x, index)
  if (is.null(index)) {
    index <- 1:(length(x) * nrow(x[[1]]))
  }
  # potential max
  doses2 <- - y[, "b2"] / (2 * y[, "b3"])
  dat1 <- posterior(
    x,
    doses = doses1,
    times = time,
    iter = index,
    return_samples = TRUE,
    return_stats = FALSE
  )$samps
  vals2 <- dreamer_mean(x, dose = doses2, time = time, index = index)
  dat2 <- data.frame(
    dose = doses2,
    iter = index,
    mean_response = vals2
  ) %>%
    dplyr::filter(.data$dose >= !!lower, .data$dose <= !!upper)
  responses <- dplyr::bind_rows(dat1, dat2)
  er <- pextreme(
    responses = responses,
    greater = greater
  )
  return(er)
}

#' @keywords internal
#' @export
get_extreme.dreamer_mcmc_quad_binary <- get_extreme.dreamer_mcmc_quad # nolint

#' @keywords internal
#' @export
get_extreme.dreamer_mcmc_logquad <- function( #nolint
  x,
  time,
  greater,
  lower,
  upper,
  index
) {
  check_time(time)
  y <- subset_mcmc(x, index)
  if (is.null(index)) {
    index <- 1:(length(x) * nrow(x[[1]]))
  }
  doses1 <- c(lower, upper)
  # possible maximum
  doses2 <- exp(- y[, "b2"] / (2 * y[, "b3"])) - 1
  dat1 <- posterior(
    x,
    doses = doses1,
    times = time,
    iter = index,
    return_samples = TRUE,
    return_stats = FALSE
  )$samps
  vals2 <- dreamer_mean(x, dose = doses2, time = time, index = index)
  dat2 <- data.frame(
    dose = doses2,
    iter = index,
    mean_response = vals2
  ) %>%
    dplyr::filter(.data$dose >= !!lower, .data$dose <= !!upper)
  responses <- dplyr::bind_rows(dat1, dat2)
  er <- pextreme(
    responses = responses,
    greater = greater
  )
  return(er)
}

#' @keywords internal
#' @export
get_extreme.dreamer_mcmc_logquad_binary <- # nolint
  get_extreme.dreamer_mcmc_logquad

#' @keywords internal
#' @export
get_extreme.dreamer_mcmc_beta <- function(
  x,
  time,
  greater,
  lower,
  upper,
  index
) {
  check_time(time)
  y <- subset_mcmc(x, index)
  if (is.null(index)) {
    index <- 1:(length(x) * nrow(x[[1]]))
  }
  scale <- attr(x, "scale")
  doses1 <- c(lower, upper)
  doses2 <- scale * y[, "b3"] / (y[, "b3"] + y[, "b4"])
  dat1 <- posterior(
    x,
    doses = doses1,
    times = time,
    iter = index,
    return_samples = TRUE,
    return_stats = FALSE
  )$samps
  vals2 <- dreamer_mean(x, dose = doses2, time = time, index = index)
  dat2 <- data.frame(
    dose = doses2,
    iter = index,
    mean_response = vals2
  ) %>%
    dplyr::filter(.data$dose >= !!lower, .data$dose <= !!upper)
  responses <- dplyr::bind_rows(dat1, dat2)
  er <- pextreme(
    responses = responses,
    greater = greater
  )
  return(er)
}

#' @keywords internal
#' @export
get_extreme.dreamer_mcmc_beta_binary <- # nolint
  get_extreme.dreamer_mcmc_beta

pextreme <- function(responses, greater) {
  dat <- responses %>%
    dplyr::group_by(.data$iter)
  if (greater) {
    out <- dat %>%
      dplyr::summarize(
        dose_ind = which.max(.data$mean_response),
        dose = .data$dose[.data$dose_ind],
        mean_response = .data$mean_response[.data$dose_ind]
      )
  } else {
    out <- dat %>%
      dplyr::summarize(
        dose_ind = which.min(.data$mean_response),
        dose = .data$dose[.data$dose_ind],
        mean_response = .data$mean_response[.data$dose_ind]
      )
  }
  out <- out %>%
    dplyr::ungroup() %>%
    dplyr::mutate(greater = !!greater) %>%
    dplyr::select(
      doses = .data$dose,
      extreme_responses = .data$mean_response,
      greater
    )
  return(out)
}

check_time <- function(time) {
  if (length(time) > 1) {
    rlang::abort("time must have length 1 or be NULL", class = "dreamer")
  }
}
