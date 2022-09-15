#' @title Generate Data from Dose Response Models
#' @description See the model definitions below for specifics for each model.
#' @return A dataframe of random subjects from the specified model and
#'   parameters.
#' @name dreamer_data
#' @param b1,b2,b3,b4 parameters in the models.  See sections below for
#'   each parameter's interpretation in a given model.
#' @param n_cohorts a vector listing the size of each cohort.
#' @param doses a vector listing the dose for each cohort.
#' @param sigma standard deviation.
#' @param link character vector indicating the link function for binary models.
#' @param times the times at which data should be simulated if a longitudinal
#'   model is specified.
#' @param t_max the t_max parameter used in the longitudinal model.
#' @param longitudinal a string indicating the longitudinal model to be used.
#'   Can be "linear", "itp", or "idp".
#' @param ... additional longitudinal parameters.
#' @inheritSection model_linear Linear
#' @inheritSection model_quad Quadratic
#' @inheritSection model_loglinear Log-linear
#' @inheritSection model_logquad Log-quadratic
#' @inheritSection model_emax EMAX
#' @inheritSection model_exp Exponential
#' @inheritSection model_linear_binary Linear Binary
#' @inheritSection model_quad_binary Quadratic Binary
#' @inheritSection model_loglinear_binary Log-linear Binary
#' @inheritSection model_logquad_binary Log-quadratic Binary
#' @inheritSection model_emax_binary EMAX Binary
#' @inheritSection model_exp_binary Exponential Binary
#' @inheritSection model_independent Independent
#' @inheritSection model_independent_binary Independent Binary
#' @inheritSection model_longitudinal Longitudinal Linear
#' @inheritSection model_longitudinal Longitudinal ITP
#' @inheritSection model_longitudinal Longitudinal IDP
#' @example man/examples/ex-dreamer_data.R
NULL

#' @describeIn dreamer_data generate data from linear dose response.
#' @export
dreamer_data_linear <- function(
  n_cohorts,
  doses,
  b1,
  b2,
  sigma,
  times,
  t_max,
  longitudinal = NULL,
  ...
) {
  check_len(doses, n_cohorts, "length(n_cohorts) must equal length(doses)")
  dat <- get_base_data(n_cohorts, doses)
  if (is.null(longitudinal)) {
    dat$response <- dreamer_mean_linear(dose = dat$dose, b1 = b1, b2 = b2)
  } else {
    dat <- tidyr::expand_grid(
      dat,
      time = times
    )
    args <- list(...)
    a <- args$a
    mcmc <- get_args(...)
    dat$response <- a +
      dreamer_mean_linear(dose = dat$dose, b1 = b1, b2 = b2) *
      dreamer_mean_longitudinal_mcmc(
        time = times,
        longitudinal_model = paste0("dreamer_longitudinal_", longitudinal),
        mcmc = mcmc,
        t_max = t_max
      )
  }
  dat$response <- dat$response + rnorm(nrow(dat), 0, sigma)
  return(dat)
}

#' @describeIn dreamer_data generate data from linear binary dose response.
#' @export
dreamer_data_linear_binary <- function(
  n_cohorts,
  doses,
  b1,
  b2,
  link,
  times,
  t_max,
  longitudinal = NULL,
  ...
) {
  ilink <- get(paste0("i", link))
  dat <- dreamer_data_linear(
    n_cohorts = n_cohorts,
    doses = doses,
    b1 = b1,
    b2 = b2,
    sigma = 0,
    times = times,
    t_max = t_max,
    longitudinal = longitudinal,
    ...
  )
  dat$response <- rbinom(nrow(dat), 1, ilink(dat$response))
  return(dat)
}

#' @describeIn dreamer_data generate data from quadratic dose response.
#' @export
dreamer_data_quad <- function(
  n_cohorts,
  doses,
  b1,
  b2,
  b3,
  sigma,
  times,
  t_max,
  longitudinal = NULL,
  ...
) {
  check_len(doses, n_cohorts, "length(n_cohorts) must equal length(doses)")
  dat <- get_base_data(n_cohorts, doses)
  if (is.null(longitudinal)) {
    dat$response <- dreamer_mean_quad(
      dose = dat$dose,
      b1 = b1,
      b2 = b2,
      b3 = b3
    )
  } else {
    dat <- tidyr::expand_grid(
      dat,
      time = times
    )
    args <- list(...)
    a <- args$a
    mcmc <- get_args(...)
    dat$response <- a +
      dreamer_mean_quad(dose = dat$dose, b1 = b1, b2 = b2, b3 = b3) *
      dreamer_mean_longitudinal_mcmc(
        time = times,
        longitudinal_model = paste0("dreamer_longitudinal_", longitudinal),
        mcmc = mcmc,
        t_max = t_max
      )
  }
  dat$response <- dat$response + rnorm(nrow(dat), 0, sigma)
  return(dat)
}

#' @describeIn dreamer_data generate data from quadratic binary dose response.
#' @export
dreamer_data_quad_binary <- function(
  n_cohorts,
  doses,
  b1,
  b2,
  b3,
  link,
  times,
  t_max,
  longitudinal = NULL,
  ...
) {
  ilink <- get(paste0("i", link))
  dat <- dreamer_data_quad(
    n_cohorts = n_cohorts,
    doses = doses,
    b1 = b1,
    b2 = b2,
    b3 = b3,
    sigma = 0,
    times = times,
    t_max = t_max,
    longitudinal = longitudinal,
    ...
  )
  dat$response <- rbinom(nrow(dat), 1, ilink(dat$response))
  return(dat)
}

#' @describeIn dreamer_data generate data from log-linear dose response.
#' @export
dreamer_data_loglinear <- function(
  n_cohorts,
  doses,
  b1,
  b2,
  sigma,
  times,
  t_max,
  longitudinal = NULL,
  ...
) {
  check_len(doses, n_cohorts, "length(n_cohorts) must equal length(doses)")
  dat <- get_base_data(n_cohorts, doses)
  if (is.null(longitudinal)) {
    dat$response <- dreamer_mean_loglinear(dose = dat$dose, b1 = b1, b2 = b2)
  } else {
    dat <- tidyr::expand_grid(
      dat,
      time = times
    )
    args <- list(...)
    a <- args$a
    mcmc <- get_args(...)
    dat$response <- a +
      dreamer_mean_loglinear(dose = dat$dose, b1 = b1, b2 = b2) *
      dreamer_mean_longitudinal_mcmc(
        time = times,
        longitudinal_model = paste0("dreamer_longitudinal_", longitudinal),
        mcmc = mcmc,
        t_max = t_max
      )
  }
  dat$response <- dat$response + rnorm(nrow(dat), 0, sigma)
  return(dat)
}

#' @describeIn dreamer_data generate data from binary log-linear dose response.
#' @export
dreamer_data_loglinear_binary <- function(
  n_cohorts,
  doses,
  b1,
  b2,
  link,
  times,
  t_max,
  longitudinal = NULL,
  ...
) {
  ilink <- get(paste0("i", link))
  dat <- dreamer_data_loglinear(
    n_cohorts = n_cohorts,
    doses = doses,
    b1 = b1,
    b2 = b2,
    sigma = 0,
    times = times,
    t_max = t_max,
    longitudinal = longitudinal,
    ...
  )
  dat$response <- rbinom(nrow(dat), 1, ilink(dat$response))
  return(dat)
}

#' @describeIn dreamer_data generate data from log-quadratic dose response.
#' @export
dreamer_data_logquad <- function(
  n_cohorts,
  doses,
  b1,
  b2,
  b3,
  sigma,
  times,
  t_max,
  longitudinal = NULL,
  ...
) {
  check_len(doses, n_cohorts, "length(n_cohorts) must equal length(doses)")
  dat <- get_base_data(n_cohorts, doses)
  if (is.null(longitudinal)) {
    dat$response <- dreamer_mean_logquad(
      dose = dat$dose,
      b1 = b1,
      b2 = b2,
      b3 = b3
    )
  } else {
    dat <- tidyr::expand_grid(
      dat,
      time = times
    )
    args <- list(...)
    a <- args$a
    mcmc <- get_args(...)
    dat$response <- a +
      dreamer_mean_logquad(dose = dat$dose, b1 = b1, b2 = b2, b3 = b3) *
      dreamer_mean_longitudinal_mcmc(
        time = times,
        longitudinal_model = paste0("dreamer_longitudinal_", longitudinal),
        mcmc = mcmc,
        t_max = t_max
      )
  }
  dat$response <- dat$response + rnorm(nrow(dat), 0, sigma)
  return(dat)
}

#' @describeIn dreamer_data generate data from log-quadratic binary dose
#'   response.
#' @export
dreamer_data_logquad_binary <- function(
  n_cohorts,
  doses,
  b1,
  b2,
  b3,
  link,
  times,
  t_max,
  longitudinal = NULL,
  ...
) {
  ilink <- get(paste0("i", link))
  dat <- dreamer_data_logquad(
    n_cohorts = n_cohorts,
    doses = doses,
    b1 = b1,
    b2 = b2,
    b3 = b3,
    sigma = 0,
    times = times,
    t_max = t_max,
    longitudinal = longitudinal,
    ...
  )
  dat$response <- rbinom(nrow(dat), 1, ilink(dat$response))
  return(dat)
}

#' @describeIn dreamer_data generate data from EMAX dose response.
#' @export
dreamer_data_emax <- function(
  n_cohorts,
  doses,
  b1,
  b2,
  b3,
  b4,
  sigma,
  times,
  t_max,
  longitudinal = NULL,
  ...
) {
  check_len(doses, n_cohorts, "length(n_cohorts) must equal length(doses)")
  dat <- get_base_data(n_cohorts, doses)
  if (is.null(longitudinal)) {
    dat$response <- dreamer_mean_emax(
      dose = dat$dose,
      b1 = b1,
      b2 = b2,
      b3 = b3,
      b4 = b4
    )
  } else {
    dat <- tidyr::expand_grid(
      dat,
      time = times
    )
    args <- list(...)
    a <- args$a
    mcmc <- get_args(...)
    dat$response <- a +
      dreamer_mean_emax(dose = dat$dose, b1 = b1, b2 = b2, b3 = b3, b4 = b4) *
      dreamer_mean_longitudinal_mcmc(
        time = times,
        longitudinal_model = paste0("dreamer_longitudinal_", longitudinal),
        mcmc = mcmc,
        t_max = t_max
      )
  }
  dat$response <- dat$response + rnorm(nrow(dat), 0, sigma)
  return(dat)
}

#' @describeIn dreamer_data generate data from EMAX binary dose response.
#' @export
dreamer_data_emax_binary <- function(
  n_cohorts,
  doses,
  b1,
  b2,
  b3,
  b4,
  link,
  times,
  t_max,
  longitudinal = NULL,
  ...
) {
  ilink <- get(paste0("i", link))
  dat <- dreamer_data_emax(
    n_cohorts = n_cohorts,
    doses = doses,
    b1 = b1,
    b2 = b2,
    b3 = b3,
    b4 = b4,
    sigma = 0,
    times = times,
    t_max = t_max,
    longitudinal = longitudinal,
    ...
  )
  dat$response <- rbinom(nrow(dat), 1, ilink(dat$response))
  return(dat)
}

#' @describeIn dreamer_data generate data from exponential dose response.
#' @export
dreamer_data_exp <- function(
  n_cohorts,
  doses,
  b1,
  b2,
  b3,
  sigma,
  times,
  t_max,
  longitudinal = NULL,
  ...
) {
  check_len(doses, n_cohorts, "length(n_cohorts) must equal length(doses)")
  dat <- get_base_data(n_cohorts, doses)
  if (is.null(longitudinal)) {
    dat$response <- dreamer_mean_exp(dose = dat$dose, b1 = b1, b2 = b2, b3 = b3)
  } else {
    dat <- tidyr::expand_grid(
      dat,
      time = times
    )
    args <- list(...)
    a <- args$a
    mcmc <- get_args(...)
    dat$response <- a +
      dreamer_mean_exp(dose = dat$dose, b1 = b1, b2 = b2, b3 = b3) *
      dreamer_mean_longitudinal_mcmc(
        time = times,
        longitudinal_model = paste0("dreamer_longitudinal_", longitudinal),
        mcmc = mcmc,
        t_max = t_max
      )
  }
  dat$response <- dat$response + rnorm(nrow(dat), 0, sigma)
  return(dat)
}

#' @describeIn dreamer_data generate data from exponential binary dose response.
#' @export
dreamer_data_exp_binary <- function(
  n_cohorts,
  doses,
  b1,
  b2,
  b3,
  link,
  times,
  t_max,
  longitudinal = NULL,
  ...
) {
  ilink <- get(paste0("i", link))
  dat <- dreamer_data_exp(
    n_cohorts = n_cohorts,
    doses = doses,
    b1 = b1,
    b2 = b2,
    b3 = b3,
    sigma = 0,
    times = times,
    t_max = t_max,
    longitudinal = longitudinal,
    ...
  )
  dat$response <- rbinom(nrow(dat), 1, ilink(dat$response))
  return(dat)
}

#' @describeIn dreamer_data generate data from Beta dose response.
#' @param scale a scaling parameter (fixed, specified by the user) for the
#'   beta models.
#' @export
dreamer_data_beta <- function(
  n_cohorts,
  doses,
  b1,
  b2,
  b3,
  b4,
  scale,
  sigma,
  times,
  t_max,
  longitudinal = NULL,
  ...
) {
  check_len(doses, n_cohorts, "length(n_cohorts) must equal length(doses)")
  dat <- get_base_data(n_cohorts, doses)
  if (is.null(longitudinal)) {
    dat$response <- dreamer_mean_beta(
      dose = dat$dose,
      b1 = b1,
      b2 = b2,
      b3 = b3,
      b4 = b4,
      scale = scale
    )
  } else {
    dat <- tidyr::expand_grid(
      dat,
      time = times
    )
    args <- list(...)
    a <- args$a
    mcmc <- get_args(...)
    dat$response <- a +
      dreamer_mean_beta(
        dose = dat$dose,
        b1 = b1,
        b2 = b2,
        b3 = b3,
        b4 = b4,
        scale = scale
      ) *
      dreamer_mean_longitudinal_mcmc(
        time = times,
        longitudinal_model = paste0("dreamer_longitudinal_", longitudinal),
        mcmc = mcmc,
        t_max = t_max
      )
  }
  dat$response <- dat$response + rnorm(nrow(dat), 0, sigma)
  return(dat)
}

#' @describeIn dreamer_data generate data from binary Beta dose response.
#' @export
dreamer_data_beta_binary <- function(
  n_cohorts,
  doses,
  b1,
  b2,
  b3,
  b4,
  scale,
  link,
  times,
  t_max,
  longitudinal = NULL,
  ...
) {
  ilink <- get(paste0("i", link))
  dat <- dreamer_data_beta(
    n_cohorts = n_cohorts,
    doses = doses,
    b1 = b1,
    b2 = b2,
    b3 = b3,
    b4 = b4,
    scale = scale,
    sigma = 0,
    times = times,
    t_max = t_max,
    longitudinal = longitudinal,
    ...
  )
  dat$response <- rbinom(nrow(dat), 1, ilink(dat$response))
  return(dat)
}

#' @describeIn dreamer_data generate data from an independent dose response.
#' @export
dreamer_data_independent <- function(
  n_cohorts,
  doses,
  b1,
  sigma,
  times,
  t_max,
  longitudinal = NULL,
  ...
) {
  check_len(doses, n_cohorts, "length(n_cohorts) must equal length(doses)")
  check_len(doses, b1, "length(b1) must equal length(doses)")
  check_duplicate_dose(doses)

  dat <- get_base_data(n_cohorts, doses)
  if (is.null(longitudinal)) {
    dose_index <- vapply(dat$dose, function(x) which(x == doses), numeric(1))
    dat$response <- dreamer_mean_independent(dose_index = dose_index, b1 = b1)
  } else {
    dat <- tidyr::expand_grid(
      dat,
      time = times
    )
    dose_index <- vapply(dat$dose, function(x) which(x == doses), numeric(1))
    args <- list(...)
    a <- args$a
    mcmc <- get_args(...)
    dat$response <- a +
      dreamer_mean_independent(dose_index = dose_index, b1 = b1) *
      dreamer_mean_longitudinal_mcmc(
        time = times,
        longitudinal_model = paste0("dreamer_longitudinal_", longitudinal),
        mcmc = mcmc,
        t_max = t_max
      )
  }
  dat$response <- dat$response + rnorm(nrow(dat), 0, sigma)
  return(dat)
}

#' @describeIn dreamer_data generate data from an independent dose response.
#' @export
dreamer_data_independent_binary <- function( #nolint
  n_cohorts,
  doses,
  b1,
  link,
  times,
  t_max,
  longitudinal = NULL,
  ...
) {
  ilink <- get(paste0("i", link))
  dat <- dreamer_data_independent(
    n_cohorts = n_cohorts,
    doses = doses,
    b1 = b1,
    sigma = 0,
    times = times,
    t_max = t_max,
    longitudinal = longitudinal,
    ...
  )
  dat$response <- rbinom(nrow(dat), 1, ilink(dat$response))
  return(dat)
}

check_len <- function(x, y, msg) {
  if (length(x) != length(y)) {
    stop(msg, call. = FALSE)
  }
}

check_duplicate_dose <- function(doses) {
  if (anyDuplicated(doses) > 0) {
    stop("duplicate doses not allowed.", call. = FALSE)
  }
}

get_args <- function(...) {
  out <- NULL
  args <- list(...)
  args$a <- NULL
  if (length(args) > 0) {
    out <- t(as.matrix(unlist(args)))
  }
  return(out)
}

get_base_data <- function(n_cohorts, doses) {
  data.frame(
    cohort = purrr::map2(
      n_cohorts,
      seq_len(length(n_cohorts)),
      ~ rep(.y, .x)
    ) %>%
      do_call("c"),
    dose = purrr::map2(n_cohorts, doses, ~ rep(.y, .x)) %>%
      do_call("c")
  ) %>%
    dplyr::mutate(subject = seq_len(n()))
}
