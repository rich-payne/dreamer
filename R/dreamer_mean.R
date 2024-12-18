dreamer_mean_linear <- function(dose, b1, b2) {
  b1 + b2 * dose
}

dreamer_mean_quad <- function(dose, b1, b2, b3) {
  b1 + b2 * dose + b3 * dose^2
}

dreamer_mean_loglinear <- function(dose, b1, b2) {
  b1 + b2 * log(dose + 1)
}

dreamer_mean_logquad <- function(dose, b1, b2, b3) {
  b1 + b2 * log(dose + 1) + b3 * log(dose + 1)^2
}

dreamer_mean_emax <- function(dose, b1, b2, b3, b4) {
  b1 + (b2 - b1) / (exp(b4 * (b3 - log(dose))) + 1)
}

dreamer_mean_exp <- function(dose, b1, b2, b3) {
  b1 + b2 * (1 - exp(- b3 * dose))
}

dreamer_mean_beta <- function(dose, b1, b2, b3, b4, scale) {
  val <- (b3 + b4) ^ (b3 + b4) / (b3 ^ b3 * b4 ^ b4)
  b1 + b2 * val * (dose / scale) ^ b3 * (1 - dose / scale) ^ b4
}

dreamer_mean_independent <- function(dose_index, b1) {
  b1[dose_index]
}

subset_mcmc <- function(x, index) {
  y <- as.matrix(x)
  if (!is.null(index)) {
    y <- y[index, , drop = FALSE]
  }
  return(y)
}

get_y <- function(y, x, index) {
  if (is.null(y)) {
    y <- subset_mcmc(x, index)
  }
  y
}

dreamer_mean <- function(x, y, dose, time, index) {
  UseMethod("dreamer_mean", x)
}

#' @export
dreamer_mean.dreamer_mcmc_linear <- function( #nolint
  x,
  y = NULL,
  dose,
  time,
  index = NULL
) {
  y <- get_y(y, x, index)
  out <- dreamer_mean_linear(dose, b1 = y[, "b1"], b2 = y[, "b2"])
  if (!is.null(attr(x, "longitudinal_model"))) {
    out <- y[, "a"] + out * dreamer_mean_longitudinal_mcmc(
      time = time,
      longitudinal_model = attr(x, "longitudinal_model"),
      mcmc = y,
      t_max = attr(x, "t_max")
    )
  }
  return(out)
}

#' @export
dreamer_mean.dreamer_mcmc_loglinear <- function( #nolint
  x,
  y = NULL,
  dose,
  time,
  index
) {
  y <- get_y(y, x, index)
  out <- dreamer_mean_loglinear(dose, b1 = y[, "b1"], b2 = y[, "b2"])
  if (!is.null(attr(x, "longitudinal_model"))) {
    out <- y[, "a"] + out * dreamer_mean_longitudinal_mcmc(
      time = time,
      longitudinal_model = attr(x, "longitudinal_model"),
      mcmc = y,
      t_max = attr(x, "t_max")
    )
  }
  return(out)
}

#' @export
dreamer_mean.dreamer_mcmc_quad <- function(x, y = NULL, dose, time, index) {
  y <- get_y(y, x, index)
  out <- dreamer_mean_quad(dose, b1 = y[, "b1"], b2 = y[, "b2"], b3 = y[, "b3"])
  if (!is.null(attr(x, "longitudinal_model"))) {
    out <- y[, "a"] + out * dreamer_mean_longitudinal_mcmc(
      time = time,
      longitudinal_model = attr(x, "longitudinal_model"),
      mcmc = y,
      t_max = attr(x, "t_max")
    )
  }
  return(out)
}

#' @export
dreamer_mean.dreamer_mcmc_logquad <- function( #nolint
  x,
  y = NULL,
  dose,
  time,
  index
) {
  y <- get_y(y, x, index)
  out <- dreamer_mean_logquad(
    dose,
    b1 = y[, "b1"],
    b2 = y[, "b2"],
    b3 = y[, "b3"]
  )
  if (!is.null(attr(x, "longitudinal_model"))) {
    out <- y[, "a"] + out * dreamer_mean_longitudinal_mcmc(
      time = time,
      longitudinal_model = attr(x, "longitudinal_model"),
      mcmc = y,
      t_max = attr(x, "t_max")
    )
  }
  return(out)
}

#' @export
dreamer_mean.dreamer_mcmc_emax <- function(x, y = NULL, dose, time, index) {
  y <- get_y(y, x, index)
  out <- dreamer_mean_emax(
    dose,
    b1 = y[, "b1"],
    b2 = y[, "b2"],
    b3 = y[, "b3"],
    b4 = y[, "b4"]
  )
  if (!is.null(attr(x, "longitudinal_model"))) {
    out <- y[, "a"] + out * dreamer_mean_longitudinal_mcmc(
      time = time,
      longitudinal_model = attr(x, "longitudinal_model"),
      mcmc = y,
      t_max = attr(x, "t_max")
    )
  }
  return(out)
}

#' @export
dreamer_mean.dreamer_mcmc_exp <- function(x, y = NULL, dose, time, index) {
  y <- get_y(y, x, index)
  out <- dreamer_mean_exp(dose, b1 = y[, "b1"], b2 = y[, "b2"], b3 = y[, "b3"])
  if (!is.null(attr(x, "longitudinal_model"))) {
    out <- y[, "a"] + out * dreamer_mean_longitudinal_mcmc(
      time = time,
      longitudinal_model = attr(x, "longitudinal_model"),
      mcmc = y,
      t_max = attr(x, "t_max")
    )
  }
  return(out)
}

#' @export
dreamer_mean.dreamer_mcmc_beta <- function(x, y = NULL, dose, time, index) {
  y <- get_y(y, x, index)
  out <- dreamer_mean_beta(
    dose,
    b1 = y[, "b1"],
    b2 = y[, "b2"],
    b3 = y[, "b3"],
    b4 = y[, "b4"],
    scale = attr(x, "scale")
  )
  if (!is.null(attr(x, "longitudinal_model"))) {
    out <- y[, "a"] +
      out * dreamer_mean_longitudinal_mcmc(
        time = time,
        longitudinal_model = attr(x, "longitudinal_model"),
        mcmc = y,
        t_max = attr(x, "t_max")
      )
  }
  return(out)
}

#' @export
dreamer_mean.dreamer_mcmc_linear_binary <- function( #nolint
  x,
  y = NULL,
  dose,
  time,
  index
) {
  link <- attr(x, "link")
  check_link(link)
  lfun <- get(paste0("i", link))
  dreamer_mean.dreamer_mcmc_linear(x, y, dose, time, index) %>%
    lfun()
}

#' @export
dreamer_mean.dreamer_mcmc_loglinear_binary <- function( #nolint
  x,
  y = NULL,
  dose,
  time,
  index
) {
  link <- attr(x, "link")
  check_link(link)
  lfun <- get(paste0("i", link))
  dreamer_mean.dreamer_mcmc_loglinear(x, y, dose, time, index) %>%
    lfun()
}

#' @export
dreamer_mean.dreamer_mcmc_quad_binary <- function( #nolint
  x,
  y = NULL,
  dose,
  time,
  index
) {
  link <- attr(x, "link")
  check_link(link)
  lfun <- get(paste0("i", link))
  dreamer_mean.dreamer_mcmc_quad(x, y, dose, time, index) %>%
    lfun()
}

#' @export
dreamer_mean.dreamer_mcmc_logquad_binary <- function( #nolint
  x,
  y = NULL,
  dose,
  time,
  index
) {
  link <- attr(x, "link")
  check_link(link)
  lfun <- get(paste0("i", link))
  dreamer_mean.dreamer_mcmc_logquad(x, y, dose, time, index) %>%
    lfun()
}

#' @export
dreamer_mean.dreamer_mcmc_emax_binary <- function( #nolint
  x,
  y = NULL,
  dose,
  time,
  index
) {
  link <- attr(x, "link")
  check_link(link)
  lfun <- get(paste0("i", link))
  dreamer_mean.dreamer_mcmc_emax(x, y, dose, time, index) %>%
    lfun()
}

#' @export
dreamer_mean.dreamer_mcmc_exp_binary <- function( #nolint
  x,
  y = NULL,
  dose,
  time,
  index
) {
  link <- attr(x, "link")
  check_link(link)
  lfun <- get(paste0("i", link))
  dreamer_mean.dreamer_mcmc_exp(x, y, dose, time, index) %>%
    lfun()
}

#' @export
dreamer_mean.dreamer_mcmc_beta_binary <- function( #nolint
  x,
  y = NULL,
  dose,
  time,
  index
) {
  link <- attr(x, "link")
  check_link(link)
  lfun <- get(paste0("i", link))
  dreamer_mean.dreamer_mcmc_beta(x, y, dose, time, index) %>%
    lfun()
}

#' @export
dreamer_mean.dreamer_mcmc_independent <- function( #nolint
  x,
  y = NULL,
  dose,
  time,
  index
) {
  y <- get_y(y, x, index)
  check_dose_length(dose)
  dose_index <- which(dose == attr(x, "doses"))
  check_dose_index(dose_index, dose, x)
  col_index <- grep("^b1", colnames(y))
  out <- y[, col_index[dose_index]]
  if (!is.null(attr(x, "longitudinal_model"))) {
    out <- y[, "a"] + out * dreamer_mean_longitudinal_mcmc(
      time = time,
      longitudinal_model = attr(x, "longitudinal_model"),
      mcmc = y,
      t_max = attr(x, "t_max")
    )
  }
  return(out)
}

#' @export
dreamer_mean.dreamer_mcmc_independent_binary <- function( #nolint
  x,
  y = NULL,
  dose,
  time,
  index
) {
  link <- attr(x, "link")
  ilink_fun <- get(paste0("i", link))
  dreamer_mean.dreamer_mcmc_independent(x, y, dose, time, index) %>%
    ilink_fun()
}

dreamer_mean_longitudinal_mcmc <- function(
  time,
  longitudinal_model,
  mcmc,
  t_max
) {
  if ("dreamer_longitudinal_linear" %in% longitudinal_model) {
    out <- time / t_max
  } else if ("dreamer_longitudinal_itp" %in% longitudinal_model) {
    out <- (1 - exp(- mcmc[, "c1"] * time)) / (1 - exp(- mcmc[, "c1"] * t_max))
  } else if ("dreamer_longitudinal_idp" %in% longitudinal_model) {
    out <-
      ((1 - exp(- mcmc[, "c1"] * time)) /
        (1 - exp(- mcmc[, "c1"] * mcmc[, "d1"]))) * (time < mcmc[, "d1"]) +
      (1 - mcmc[, "gam"] * ((1 - exp(- mcmc[, "c2"] * (time - mcmc[, "d1"]))) /
        (1 - exp(- mcmc[, "c2"] * (mcmc[, "d2"] - mcmc[, "d1"]))))) *
        ((time >= mcmc[, "d1"]) & (time <= mcmc[, "d2"])) +
      (1 - mcmc[, "gam"]) * (time > mcmc[, "d2"])
  } else if (is.null(longitudinal_model)) {
    out <- 1
  } else {
    stop(paste0("Longitudinal model '", longitudinal_model, "' not supported."))
  }
  return(out)
}

check_link <- function(link) {
  if (!(link %in% c("probit", "logit"))) {
    stop("link must be 'probit' or 'logit'")
  }
}

check_dose_length <- function(dose) {
  if (length(dose) != 1) {
    stop("dose must have length 1.", call. = FALSE)
  }
}

check_dose_index <- function(dose_index, dose, x) {
  if (length(dose_index) == 0) {
    stop(
      paste0(
        "Must choose a dose (you chose ", dose, ") from existing doses: ",
        paste(attr(x, "doses"), collapse = ", ")
      ),
      call. = FALSE
    )
  }
}
