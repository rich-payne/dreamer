get_dose <- function(x, time, response, index, lower, upper) {
  UseMethod("get_dose", x)
}

# get the lowest dose that have a certain effect within interval [lower, upper]

#' @keywords internal
#' @export
get_dose.dreamer_mcmc_linear <- function(
  x,
  time,
  response,
  index,
  lower,
  upper
) {
  check_time(time)
  y <- get_data(x, index)
  mean_time <- dreamer_mean_longitudinal_mcmc(
    time = time,
    longitudinal_model = attr(x, "longitudinal_model"),
    mcmc = y,
    t_max = attr(x, "t_max")
  )
  a <- set_a(y)
  doses <- (response - y[, "b1"] * mean_time - a) /
    (y[, "b2"] * mean_time)
  doses[doses < lower | doses > upper] <- NA
  return(doses)
}

#' @keywords internal
#' @export
get_dose.dreamer_mcmc_loglinear <- function( #nolint
  x,
  time,
  response,
  index,
  lower,
  upper
) {
  check_time(time)
  y <- get_data(x, index)
  mean_time <- dreamer_mean_longitudinal_mcmc(
    time = time,
    longitudinal_model = attr(x, "longitudinal_model"),
    mcmc = y,
    t_max = attr(x, "t_max")
  )
  a <- set_a(y)
  doses <- exp(
    (response - y[, "b1"] * mean_time - a) / (y[, "b2"] * mean_time)
  ) - 1
  doses[doses < lower | doses > upper] <- NA
  return(doses)
}

#' @keywords internal
#' @export
get_dose.dreamer_mcmc_quad <- function(
  x,
  time,
  response,
  index,
  lower,
  upper
) {
  check_time(time)
  y <- get_data(x, index)
  mean_time <- dreamer_mean_longitudinal_mcmc(
    time = time,
    longitudinal_model = attr(x, "longitudinal_model"),
    mcmc = y,
    t_max = attr(x, "t_max")
  )
  mean_time <- expand_mean_time(mean_time, y)
  a <- set_a(y)
  doses <- vapply(
    seq_len(nrow(y)),
    function(i) {
      solve_quad(
        y[i, "b3"] * mean_time[i],
        y[i, "b2"] * mean_time[i],
        y[i, "b1"] * mean_time[i] - response[i] + a[i]
      )
    },
    numeric(2)
  ) %>%
    get_min_dose(lower = lower, upper = upper)
  return(doses)
}

#' @keywords internal
#' @export
get_dose.dreamer_mcmc_logquad <- function(
  x,
  time,
  response,
  index,
  lower,
  upper
) {
  check_time(time)
  y <- get_data(x, index)
  mean_time <- dreamer_mean_longitudinal_mcmc(
    time = time,
    longitudinal_model = attr(x, "longitudinal_model"),
    mcmc = y,
    t_max = attr(x, "t_max")
  )
  mean_time <- expand_mean_time(mean_time, y)
  a <- set_a(y)
  doses <- vapply(
    seq_len(nrow(y)),
    function(i) {
      solve_quad(
        y[i, "b3"] * mean_time[i],
        y[i, "b2"] * mean_time[i],
        y[i, "b1"] * mean_time[i] - response[i] + a[i]
      )
    },
    numeric(2)
  )
  doses <- exp(doses) - 1
  doses <- get_min_dose(doses = doses, lower = lower, upper = upper)
  return(doses)
}

#' @keywords internal
#' @export
get_dose.dreamer_mcmc_emax <- function(
  x,
  time,
  response,
  index,
  lower,
  upper
) {
  check_time(time)
  y <- get_data(x, index)
  mean_time <- dreamer_mean_longitudinal_mcmc(
    time = time,
    longitudinal_model = attr(x, "longitudinal_model"),
    mcmc = y,
    t_max = attr(x, "t_max")
  )
  a <- set_a(y)
  b1 <- y[, "b1"]
  b2 <- y[, "b2"]
  b3 <- y[, "b3"]
  b4 <- y[, "b4"]
  doses <- (mean_time * (b2 - b1) /
              (response - b1 * mean_time - a) - 1) ^ (-1 / b4) * exp(b3)
  doses[doses < lower | doses > upper] <- NA
  return(doses)
}

#' @keywords internal
#' @export
get_dose.dreamer_mcmc_exp <- function(
  x,
  time,
  response,
  index,
  lower,
  upper
) {
  check_time(time)
  y <- get_data(x, index)
  mean_time <- dreamer_mean_longitudinal_mcmc(
    time = time,
    longitudinal_model = attr(x, "longitudinal_model"),
    mcmc = y,
    t_max = attr(x, "t_max")
  )
  a <- set_a(y)
  b1 <- y[, "b1"] * mean_time
  b2 <- y[, "b2"] * mean_time
  b3 <- y[, "b3"]
  z <- 1 - (response - a - b1) / b2
  ind <- z > 0
  doses_positive <- - 1 / b3[ind] * log(z[ind])
  doses_positive[doses_positive < lower | doses_positive > upper] <- NA
  doses <- rep(NA_real_, nrow(y))
  doses[ind] <- doses_positive
  return(doses)
}

#' @keywords internal
#' @export
get_dose.dreamer_mcmc_beta <- function(
  x,
  time,
  response,
  index,
  lower,
  upper
) {
  check_time(time)
  y <- get_data(x, index)
  mean_time <- dreamer_mean_longitudinal_mcmc(
    time = time,
    longitudinal_model = attr(x, "longitudinal_model"),
    mcmc = y,
    t_max = attr(x, "t_max")
  )
  mean_time <- expand_mean_time(mean_time, y)
  f_opt <- function(dose, b1, b2, b3, b4, scale, response) {
    dreamer_mean_beta(
      dose = dose,
      b1 = b1,
      b2 = b2,
      b3 = b3,
      b4 = b4,
      scale = scale
    ) - response
  }
  a <- set_a(y)
  doses <- vapply(
    seq_len(nrow(y)),
    function(ind) {
      doses <- rootSolve::uniroot.all(
        f = f_opt,
        lower = lower,
        upper = upper,
        b1 = y[ind, "b1"] * mean_time[ind],
        b2 = y[ind, "b2"] * mean_time[ind],
        b3 = y[ind, "b3"],
        b4 = y[ind, "b4"],
        scale = attr(x, "scale"),
        response = response[ind] - a[ind]
      )
      get_min_na(doses)
    },
    numeric(1)
  )
}

#' @keywords internal
#' @export
get_dose.dreamer_mcmc_linear_binary <- function( #nolint
  x,
  time,
  response,
  index,
  lower,
  upper
) {
  check_time(time)
  link <- attr(x, "link")
  lfun <- get(link)
  response <- lfun(response)
  get_dose.dreamer_mcmc_linear(
    x = x,
    time = time,
    response = response,
    index = index,
    lower = lower,
    upper = upper
  )
}

#' @keywords internal
#' @export
get_dose.dreamer_mcmc_loglinear_binary <- function( #nolint
  x,
  time,
  response,
  index,
  lower,
  upper
) {
  check_time(time)
  link <- attr(x, "link")
  lfun <- get(link)
  response <- lfun(response)
  get_dose.dreamer_mcmc_loglinear(
    x = x,
    time = time,
    response = response,
    index = index,
    lower = lower,
    upper = upper
  )
}

#' @keywords internal
#' @export
get_dose.dreamer_mcmc_quad_binary <- function( #nolint
  x,
  time,
  response,
  index,
  lower,
  upper
) {
  check_time(time)
  link <- attr(x, "link")
  lfun <- get(link)
  response <- lfun(response)
  get_dose.dreamer_mcmc_quad(
    x = x,
    time = time,
    response = response,
    index = index,
    lower = lower,
    upper = upper
  )
}

#' @keywords internal
#' @export
get_dose.dreamer_mcmc_logquad_binary <- function( #nolint
  x,
  time,
  response,
  index,
  lower,
  upper
) {
  check_time(time)
  link <- attr(x, "link")
  lfun <- get(link)
  response <- lfun(response)
  get_dose.dreamer_mcmc_logquad(
    x = x,
    time = time,
    response = response,
    index = index,
    lower = lower,
    upper = upper
  )
}

#' @keywords internal
#' @export
get_dose.dreamer_mcmc_emax_binary <- function( #nolint
  x,
  time,
  response,
  index,
  lower,
  upper
) {
  check_time(time)
  link <- attr(x, "link")
  lfun <- get(link)
  response <- lfun(response)
  get_dose.dreamer_mcmc_emax(
    x = x,
    time = time,
    response = response,
    index = index,
    lower = lower,
    upper = upper
  )
}

#' @keywords internal
#' @export
get_dose.dreamer_mcmc_exp_binary <- function( #nolint
  x,
  time,
  response,
  index,
  lower,
  upper
) {
  check_time(time)
  link <- attr(x, "link")
  lfun <- get(link)
  response <- lfun(response)
  get_dose.dreamer_mcmc_exp(
    x = x,
    time = time,
    response = response,
    index = index,
    lower = lower,
    upper = upper
  )
}

#' @keywords internal
#' @export
get_dose.dreamer_mcmc_beta_binary <- function( #nolint
  x,
  time,
  response,
  index,
  lower,
  upper
) {
  check_time(time)
  link <- attr(x, "link")
  lfun <- get(link)
  response <- lfun(response)
  get_dose.dreamer_mcmc_beta(
    x = x,
    time = time,
    response = response,
    index = index,
    lower = lower,
    upper = upper
  )
}

get_data <- function(x, index) {
  out <- as.matrix(x)[index, , drop = FALSE]
  return(out)
}

solve_quad <- function(a, b, c) {
  z <- b ^ 2 - 4 * a * c
  if (z >= 0) {
    out <- (-b + c(-1, 1) * sqrt(z)) / (2 * a)
  } else {
    out <- rep(NA_real_, 2)
  }
  return(out)
}

get_min_dose <- function(doses, lower, upper) {
  apply(
    doses,
    2,
    function(xx) {
      if (any(!is.na(xx))) {
        xx <- xx[!is.na(xx)]
        cands <- xx[xx >= lower & xx <= upper]
        return(get_min_na(cands))
      } else {
        return(NA_real_)
      }
    }
  )
}

get_min_na <- function(x) {
  if (length(x) > 0) {
    return(min(x))
  } else {
    return(NA_real_)
  }
}

expand_mean_time <- function(mean_time, mcmc) {
  if ((length(mean_time) == 1) & (nrow(mcmc) > 1)) {
    return(rep(mean_time, nrow(mcmc)))
  } else{
    return(mean_time)
  }
}

set_a <- function(y) {
  a <- rep(0, nrow(y))
  if (any(colnames(y) == "a")) {
    a <- y[, "a"]
  }
  return(a)
}
