#' @title Model Creation
#' @description Functions which set the hyperparameters, seeds, and prior
#'   weight for each model to be used in Bayesian model averaging
#'   via `dreamer_mcmc()`.
#'
#'   See each function's section below for the model's details.  In the
#'   following, \eqn{y} denotes the response variable and \eqn{d} represents
#'   the dose.
#'
#'   For the longitudinal specifications, see documentation on
#'   \code{\link{model_longitudinal}}.
#' @param w_prior a scalar between 0 and 1 indicating the prior weight of the
#'   model.
#' @param link a character string of either "logit" or "probit" indicating
#'   the link function for binary model.
#' @param longitudinal output from a call to one of the model_longitudinal_*()
#'   functions.  This is used to specify a longitudinal dose-response model.
#' @param mu_b1,sigma_b1,mu_b2,sigma_b2,mu_b3,sigma_b3,mu_b4,sigma_b4,shape,rate
#'   models parameters.  See sections below for interpretation in
#'   specific models.
#' @return A named list of the arguments in the function call.  The list has
#'   S3 classes assigned which are used internally within `dreamer_mcmc()`.
#' @name model
#' @inheritSection model_longitudinal Longitudinal Linear
#' @inheritSection model_longitudinal Longitudinal ITP
#' @inheritSection model_longitudinal Longitudinal IDP
#' @example man/examples/ex-dreamer_mcmc.R
NULL

#' @rdname model
#' @section Linear:
#'   \deqn{y \sim N(f(d), \sigma^2)}
#'   \deqn{f(d) = b_1 + b_2 * d}
#'   \deqn{b_1 \sim N(mu_b1, sigma_b1 ^ 2)}
#'   \deqn{b_2 \sim N(mu_b2, sigma_b2 ^ 2)}
#'   \deqn{1 / \sigma^2 \sim Gamma(shape, rate)}
#' @export
model_linear <- function(
  mu_b1,
  sigma_b1,
  mu_b2,
  sigma_b2,
  shape,
  rate,
  w_prior = 1,
  longitudinal = NULL
) {
  mod <- list(
    mu_b1 = mu_b1,
    sigma_b1 = sigma_b1,
    mu_b2 = mu_b2,
    sigma_b2 = sigma_b2,
    shape = shape,
    rate = rate,
    w_prior = w_prior,
    longitudinal = longitudinal
  )
  class(mod) <- c("dreamer_linear", "dreamer_continuous")
  attr(mod, "type") <- "linear"
  return(mod)
}

#' @rdname model
#' @section Quadratic:
#'   \deqn{y \sim N(f(d), \sigma^2)}
#'   \deqn{f(d) = b_1 + b_2 * d + b_3 * d^2}
#'   \deqn{b_1 \sim N(mu_b1, sigma_b1 ^ 2)}
#'   \deqn{b_2 \sim N(mu_b2, sigma_b2 ^ 2)}
#'   \deqn{b_3 \sim N(mu_b3, sigma_b3 ^ 2)}
#'   \deqn{1 / \sigma^2 \sim Gamma(shape, rate)}
#' @export
model_quad <- function(
  mu_b1,
  sigma_b1,
  mu_b2,
  sigma_b2,
  mu_b3,
  sigma_b3,
  shape,
  rate,
  w_prior = 1,
  longitudinal = NULL
) {
  mod <- list(
    mu_b1 = mu_b1,
    sigma_b1 = sigma_b1,
    mu_b2 = mu_b2,
    sigma_b2 = sigma_b2,
    mu_b3 = mu_b3,
    sigma_b3 = sigma_b3,
    shape = shape,
    rate = rate,
    w_prior = w_prior,
    longitudinal = longitudinal
  )
  class(mod) <- c("dreamer_quad", "dreamer_continuous")
  attr(mod, "type") <- "quadratic"
  return(mod)
}

#' @rdname model
#' @section Log-linear:
#'   \deqn{y \sim N(f(d), \sigma^2)}
#'   \deqn{f(d) = b_1 + b_2 * log(d + 1)}
#'   \deqn{b_1 \sim N(mu_b1, sigma_b1 ^ 2)}
#'   \deqn{b_2 \sim N(mu_b2, sigma_b2 ^ 2)}
#'   \deqn{1 / \sigma^2 \sim Gamma(shape, rate)}
#' @export
model_loglinear <- function(
  mu_b1,
  sigma_b1,
  mu_b2,
  sigma_b2,
  shape,
  rate,
  w_prior = 1,
  longitudinal = NULL
) {
  mod <- list(
    mu_b1 = mu_b1,
    sigma_b1 = sigma_b1,
    mu_b2 = mu_b2,
    sigma_b2 = sigma_b2,
    shape = shape,
    rate = rate,
    w_prior = w_prior,
    longitudinal = longitudinal
  )
  class(mod) <- c("dreamer_loglinear", "dreamer_continuous")
  attr(mod, "type") <- "log-linear"
  return(mod)
}

#' @rdname model
#' @section Log-quadratic:
#'   \deqn{y \sim N(f(d), \sigma^2)}
#'   \deqn{f(d) = b_1 + b_2 * log(d + 1) + b_3 * log(d + 1)^2}
#'   \deqn{b_1 \sim N(mu_b1, sigma_b1 ^ 2)}
#'   \deqn{b_2 \sim N(mu_b2, sigma_b2 ^ 2)}
#'   \deqn{b_3 \sim N(mu_b3, sigma_b3 ^ 2)}
#'   \deqn{1 / \sigma^2 \sim Gamma(shape, rate)}
#' @export
model_logquad <- function(
  mu_b1,
  sigma_b1,
  mu_b2,
  sigma_b2,
  mu_b3,
  sigma_b3,
  shape,
  rate,
  w_prior = 1,
  longitudinal = NULL
) {
  mod <- list(
    mu_b1 = mu_b1,
    sigma_b1 = sigma_b1,
    mu_b2 = mu_b2,
    sigma_b2 = sigma_b2,
    mu_b3 = mu_b3,
    sigma_b3 = sigma_b3,
    shape = shape,
    rate = rate,
    w_prior = w_prior,
    longitudinal = longitudinal
  )
  class(mod) <- c("dreamer_logquad", "dreamer_continuous")
  attr(mod, "type") <- "log-quadratic"
  return(mod)
}

#' @rdname model
#' @section EMAX:
#'   \deqn{y \sim N(f(d), \sigma^2)}
#'   \deqn{f(d) = b1 + b2 * d ^ b4 / (b3 ^ b4 + d ^ b4)}
#'   \deqn{b_1 \sim N(mu_b1, sigma_b1 ^ 2)}
#'   \deqn{b_2 \sim N(mu_b2, sigma_b2 ^ 2)}
#'   \deqn{b_3 \sim N(mu_b3, sigma_b3 ^ 2)}
#'   \deqn{b_4 \sim N(mu_b4, sigma_b4 ^ 2), (Truncated above 0)}
#'   \deqn{1 / \sigma^2 \sim Gamma(shape, rate)}
#'   Here, \eqn{b_1} is the placebo effect (dose = 0), \eqn{b_2} is the
#'   maximum treatment effect, \eqn{b_3} is the \eqn{ED_50}, and
#'   \eqn{b_4} is the hill or rate parameter.
#' @export
model_emax <- function(
  mu_b1,
  sigma_b1,
  mu_b2,
  sigma_b2,
  mu_b3,
  sigma_b3,
  mu_b4,
  sigma_b4,
  shape,
  rate,
  w_prior = 1,
  longitudinal = NULL
) {
  mod <- list(
    mu_b1 = mu_b1,
    sigma_b1 = sigma_b1,
    mu_b2 = mu_b2,
    sigma_b2 = sigma_b2,
    mu_b3 = mu_b3,
    sigma_b3 = sigma_b3,
    mu_b4 = mu_b4,
    sigma_b4 = sigma_b4,
    shape = shape,
    rate = rate,
    w_prior = w_prior,
    longitudinal = longitudinal
  )
  class(mod) <- c("dreamer_emax", "dreamer_continuous")
  attr(mod, "type") <- "EMAX"
  return(mod)
}

#' @rdname model
#' @param bound is the upper truncation bound of the exponential scale parameter \eqn{b_3}
#' @section Exponential:
#'   \deqn{y \sim N(f(d), \sigma^2)}
#'   \deqn{f(d) = b_1 + b_2 * (exp(d / b_3) - 1)}
#'   \deqn{b_1 \sim N(mu_b1, sigma_b1 ^ 2)}
#'   \deqn{b_2 \sim N(mu_b2, sigma_b2 ^ 2)}
#'   \deqn{b_3 \sim N(mu_b3, sigma_b3 ^ 2) T(0, bound)}
#'   \deqn{1 / \sigma^2 \sim Gamma(shape, rate)}
#' If sampling error for \eqn{b_3} then lower upper truncation bound
#' @export
model_exp <- function(
  mu_b1,
  sigma_b1,
  mu_b2,
  sigma_b2,
  mu_b3,
  sigma_b3,
  shape,
  rate,
  bound = 10,
  w_prior = 1,
  longitudinal = NULL
) {
  mod <- list(
    mu_b1 = mu_b1,
    sigma_b1 = sigma_b1,
    mu_b2 = mu_b2,
    sigma_b2 = sigma_b2,
    mu_b3 = mu_b3,
    sigma_b3 = sigma_b3,
    shape = shape,
    rate = rate,
    bound = bound,
    w_prior = w_prior,
    longitudinal = longitudinal
  )
  class(mod) <- c("dreamer_exp", "dreamer_continuous")
  attr(mod, "type") <- "exponential"
  return(mod)
}

#' @rdname model
#' @param scale a scale parameter in the Beta model. Default is 1.2 * max(dose).
#' @section Beta:
#'   \deqn{y \sim N(f(d), \sigma^2)}
#'   \deqn{f(d) = b_1 + b_2 * ((b3 + b4) ^ (b3 + b4)) /
#'     (b3 ^ b3 * b4 ^ b4) * (d / scale) ^ b3 *
#'     (1 - d / scale) ^ b4}
#'   \deqn{b_1 \sim N(mu_b1, sigma_b1 ^ 2)}
#'   \deqn{b_2 \sim N(mu_b2, sigma_b2 ^ 2)}
#'   \deqn{b_3 \sim N(mu_b3, sigma_b3 ^ 2), Truncated above 0}
#'   \deqn{b_4 \sim N(mu_b4, sigma_b4 ^ 2), Truncated above 0}
#'   \deqn{1 / \sigma^2 \sim Gamma(shape, rate)}
#'   Note that \eqn{scale} is a hyperparameter specified by the
#'   user.
#' @export
model_beta <- function(
  mu_b1,
  sigma_b1,
  mu_b2,
  sigma_b2,
  mu_b3,
  sigma_b3,
  mu_b4,
  sigma_b4,
  shape,
  rate,
  scale = NULL,
  w_prior = 1,
  longitudinal = NULL
) {
  mod <- list(
    mu_b1 = mu_b1,
    sigma_b1 = sigma_b1,
    mu_b2 = mu_b2,
    sigma_b2 = sigma_b2,
    mu_b3 = mu_b3,
    sigma_b3 = sigma_b3,
    mu_b4 = mu_b4,
    sigma_b4 = sigma_b4,
    shape = shape,
    rate = rate,
    scale = scale,
    w_prior = w_prior,
    longitudinal = longitudinal
  )
  class(mod) <- c("dreamer_beta", "dreamer_continuous")
  attr(mod, "type") <- "beta"
  return(mod)
}

#' @rdname model
#' @param doses the doses in the dataset to be modeled.  The order of the
#'   doses corresponds to the order in which the priors are specified in
#'   `mu_b1` and `sigma_b1`.
#' @section Independent:
#'   \deqn{y \sim N(f(d), \sigma^2)}
#'   \deqn{f(d) = b_{1d}}
#'   \deqn{b_{1d} \sim N(mu_b1[d], sigma_b1[d] ^ 2)}
#'   \deqn{1 / \sigma^2 \sim Gamma(shape, rate)}
#' @section Independent Details:
#'   The independent model models the effect of each dose independently.
#'   Vectors can be supplied to `mu_b1` and `sigma_b1` to set a different
#'   prior for each dose in the order the doses are supplied to `doses`.
#'   If scalars are supplied to `mu_b1` and `sigma_b1`, then the same prior
#'   is used for each dose, and the `doses` argument is not needed.
#' @export
model_independent <- function(
  mu_b1,
  sigma_b1,
  shape,
  rate,
  doses = NULL,
  w_prior = 1,
  longitudinal = NULL
) {
  check_ind_model_parms(doses, mu_b1, sigma_b1)
  mod <- list(
    mu_b1 = mu_b1,
    sigma_b1 = sigma_b1,
    shape = shape,
    rate = rate,
    doses = doses,
    w_prior = w_prior,
    longitudinal = longitudinal
  )
  class(mod) <- c(
    "dreamer_independent_continuous",
    "dreamer_independent",
    "dreamer_continuous"
  )
  attr(mod, "type") <- "independent"
  return(mod)
}

#' @rdname model
#' @section Linear Binary:
#'   \deqn{y \sim Binomial(n, f(d))}
#'   \deqn{link(f(d)) = b_1 + b_2 * d}
#'   \deqn{b_1 \sim N(mu_b1, sigma_b1 ^ 2)}
#'   \deqn{b_2 \sim N(mu_b2, sigma_b2 ^ 2)}
#' @export
model_linear_binary <- function(
  mu_b1,
  sigma_b1,
  mu_b2,
  sigma_b2,
  link,
  w_prior = 1,
  longitudinal = NULL
) {
  mod <- list(
    mu_b1 = mu_b1,
    sigma_b1 = sigma_b1,
    mu_b2 = mu_b2,
    sigma_b2 = sigma_b2,
    w_prior = w_prior,
    link = link,
    longitudinal = longitudinal
  )
  class(mod) <- c("dreamer_linear_binary", "dreamer_binary")
  attr(mod, "type") <- "linear"
  return(mod)
}

#' @rdname model
#' @section Quadratic Binary:
#'   \deqn{y \sim Binomial(n, f(d))}
#'   \deqn{link(f(d)) = b_1 + b_2 * d + b_3 * d^2}
#'   \deqn{b_1 \sim N(mu_b1, sigma_b1 ^ 2)}
#'   \deqn{b_2 \sim N(mu_b2, sigma_b2 ^ 2)}
#'   \deqn{b_3 \sim N(mu_b3, sigma_b3 ^ 2)}
#' @export
model_quad_binary <- function(
  mu_b1,
  sigma_b1,
  mu_b2,
  sigma_b2,
  mu_b3,
  sigma_b3,
  link,
  w_prior = 1,
  longitudinal = NULL
) {
  mod <- list(
    mu_b1 = mu_b1,
    sigma_b1 = sigma_b1,
    mu_b2 = mu_b2,
    sigma_b2 = sigma_b2,
    mu_b3 = mu_b3,
    sigma_b3 = sigma_b3,
    link = link,
    w_prior = w_prior,
    longitudinal = longitudinal
  )
  class(mod) <- c("dreamer_quad_binary", "dreamer_binary")
  attr(mod, "type") <- "quadratic"
  return(mod)
}

#' @rdname model
#' @section Log-linear Binary:
#'   \deqn{y \sim Binomial(n, f(d))}
#'   \deqn{link(f(d)) = b_1 + b_2 * log(d + 1)}
#'   \deqn{b_1 \sim N(mu_b1, sigma_b1 ^ 2)}
#'   \deqn{b_2 \sim N(mu_b2, sigma_b2 ^ 2)}
#' @export
model_loglinear_binary <- function(
  mu_b1,
  sigma_b1,
  mu_b2,
  sigma_b2,
  link,
  w_prior = 1,
  longitudinal = NULL
) {
  mod <- list(
    mu_b1 = mu_b1,
    sigma_b1 = sigma_b1,
    mu_b2 = mu_b2,
    sigma_b2 = sigma_b2,
    link = link,
    w_prior = w_prior,
    longitudinal = longitudinal
  )
  class(mod) <- c("dreamer_loglinear_binary", "dreamer_binary")
  attr(mod, "type") <- "log-linear"
  return(mod)
}

#' @rdname model
#' @section Log-quadratic Binary:
#'   \deqn{y \sim Binomial(n, f(d))}
#'   \deqn{link(f(d)) = b_1 + b_2 * log(d + 1) + b_3 * log(d + 1)^2}
#'   \deqn{b_1 \sim N(mu_b1, sigma_b1 ^ 2)}
#'   \deqn{b_2 \sim N(mu_b2, sigma_b2 ^ 2)}
#'   \deqn{b_3 \sim N(mu_b3, sigma_b3 ^ 2)}
#' @export
model_logquad_binary <- function(
  mu_b1,
  sigma_b1,
  mu_b2,
  sigma_b2,
  mu_b3,
  sigma_b3,
  link,
  w_prior = 1,
  longitudinal = NULL
) {
  mod <- list(
    mu_b1 = mu_b1,
    sigma_b1 = sigma_b1,
    mu_b2 = mu_b2,
    sigma_b2 = sigma_b2,
    mu_b3 = mu_b3,
    sigma_b3 = sigma_b3,
    link = link,
    w_prior = w_prior,
    longitudinal = longitudinal
  )
  class(mod) <- c("dreamer_logquad_binary", "dreamer_binary")
  attr(mod, "type") <- "log-quadratic"
  return(mod)
}

#' @rdname model
#' @section EMAX Binary:
#'   \deqn{y \sim Binomial(n, f(d))}
#'   \deqn{link(f(d)) = b1 + b2 * d ^ b4 / (b3 ^ b4 + d ^ b4)}
#'   \deqn{b_1 \sim N(mu_b1, sigma_b1 ^ 2)}
#'   \deqn{b_2 \sim N(mu_b2, sigma_b2 ^ 2)}
#'   \deqn{b_3 \sim N(mu_b3, sigma_b3 ^ 2)}
#'   \deqn{b_4 \sim N(mu_b4, sigma_b4 ^ 2), (Truncated above 0)}
#'   Here, on the \eqn{link(f(d))} scale,
#'   \eqn{b_1} is the placebo effect (dose = 0), \eqn{b_2} is the
#'   maximum treatment effect, \eqn{b_3} is the \eqn{ED_50}, and
#'   \eqn{b_4} is the hill or rate parameter.
#' @export
model_emax_binary <- function(
  mu_b1,
  sigma_b1,
  mu_b2,
  sigma_b2,
  mu_b3,
  sigma_b3,
  mu_b4,
  sigma_b4,
  link,
  w_prior = 1,
  longitudinal = NULL
) {
  mod <- list(
    mu_b1 = mu_b1,
    sigma_b1 = sigma_b1,
    mu_b2 = mu_b2,
    sigma_b2 = sigma_b2,
    mu_b3 = mu_b3,
    sigma_b3 = sigma_b3,
    mu_b4 = mu_b4,
    sigma_b4 = sigma_b4,
    link = link,
    w_prior = w_prior,
    longitudinal = longitudinal
  )
  class(mod) <- c("dreamer_emax_binary", "dreamer_binary")
  attr(mod, "type") <- "EMAX"
  return(mod)
}

#' @rdname model
#' @param bound is the upper truncation bound of the exponential scale parameter \eqn{b_3}
#' @section Exponential Binary:
#'   \deqn{y \sim N(f(d), \sigma^2)}
#'   \deqn{f(d) = b_1 + b_2 * (exp(d / b_3) - 1)}
#'   \deqn{b_1 \sim N(mu_b1, sigma_b1 ^ 2)}
#'   \deqn{b_2 \sim N(mu_b2, sigma_b2 ^ 2)}
#'   \deqn{b_3 \sim N(mu_b3, sigma_b3 ^ 2) T(0, bound)}
#' If sampling error for \eqn{b_3} then lower upper truncation bound
#' @export
model_exp_binary <- function(
  mu_b1,
  sigma_b1,
  mu_b2,
  sigma_b2,
  mu_b3,
  sigma_b3,
  bound = 10,
  link,
  w_prior = 1,
  longitudinal = NULL
) {
  mod <- list(
    mu_b1 = mu_b1,
    sigma_b1 = sigma_b1,
    mu_b2 = mu_b2,
    sigma_b2 = sigma_b2,
    mu_b3 = mu_b3,
    sigma_b3 = sigma_b3,
    bound = bound,
    link = link,
    w_prior = w_prior,
    longitudinal = longitudinal
  )
  class(mod) <- c("dreamer_exp_binary", "dreamer_binary")
  attr(mod, "type") <- "exponential"
  return(mod)
}

#' @rdname model
#' @section Beta Binary:
#'   \deqn{y \sim Binomial(n, f(d)}
#'   \deqn{link(f(d)) = b_1 + b_2 * ((b3 + b4) ^ (b3 + b4)) /
#'     (b3 ^ b3 * b4 ^ b4) * (d / scale) ^ b3 *
#'     (1 - d / scale) ^ b4}
#'   \deqn{b_1 \sim N(mu_b1, sigma_b1 ^ 2)}
#'   \deqn{b_2 \sim N(mu_b2, sigma_b2 ^ 2)}
#'   \deqn{b_3 \sim N(mu_b3, sigma_b3 ^ 2), Truncated above 0}
#'   \deqn{b_4 \sim N(mu_b4, sigma_b4 ^ 2), Truncated above 0}
#'   Note that \eqn{scale} is a hyperparameter specified by the
#'   user.
#' @export
model_beta_binary <- function(
  mu_b1,
  sigma_b1,
  mu_b2,
  sigma_b2,
  mu_b3,
  sigma_b3,
  mu_b4,
  sigma_b4,
  scale = NULL,
  link,
  w_prior = 1,
  longitudinal = NULL
) {
  mod <- list(
    mu_b1 = mu_b1,
    sigma_b1 = sigma_b1,
    mu_b2 = mu_b2,
    sigma_b2 = sigma_b2,
    mu_b3 = mu_b3,
    sigma_b3 = sigma_b3,
    mu_b4 = mu_b4,
    sigma_b4 = sigma_b4,
    scale = scale,
    link = link,
    w_prior = w_prior,
    longitudinal = longitudinal
  )
  class(mod) <- c("dreamer_beta_binary", "dreamer_binary")
  attr(mod, "type") <- "beta"
  return(mod)
}

#' @rdname model
#' @section Independent Binary:
#'   \deqn{y \sim Binomial(n, f(d))}
#'   \deqn{link(f(d)) = b_{1d}}
#'   \deqn{b_{1d} \sim N(mu_b1[d], sigma_b1[d]) ^ 2}
#' @section Independent Binary Details:
#'   The independent model models the effect of each dose independently.
#'   Vectors can be supplied to `mu_b1` and `sigma_b1` to set a different
#'   prior for each dose in the order the doses are supplied to `doses`.
#'   If scalars are supplied to `mu_b1` and `sigma_b1`, then the same prior
#'   is used for each dose, and the `doses` argument is not needed.
#' @export
model_independent_binary <- function(
  mu_b1,
  sigma_b1,
  doses = NULL,
  link,
  w_prior = 1,
  longitudinal = NULL
) {
  check_ind_model_parms(doses, mu_b1, sigma_b1)
  mod <- list(
    mu_b1 = mu_b1,
    sigma_b1 = sigma_b1,
    doses = doses,
    link = link,
    w_prior = w_prior,
    longitudinal = longitudinal
  )
  class(mod) <- c(
    "dreamer_independent_binary",
    "dreamer_independent",
    "dreamer_binary"
  )
  attr(mod, "type") <- "independent"
  return(mod)
}

check_ind_model_parms <- function(doses, mu_b1, sigma_b1) {
  if (length(mu_b1) != length(sigma_b1)) {
    stop("length(mu_b1) must equal length(sigma_b1)")
  }
  if (length(mu_b1) > 1 && is.null(doses)) {
    stop("doses must be specified if mu_b1 and sigma_b1 have length > 1")
  }
  if (length(mu_b1) > 1 && (length(mu_b1) != length(doses))) {
    stop(
      "lengths of mu_b1, sigma_b1, and doses must be equal if length(mu_b1) > 1"
    )
  }
  if (!is.null(doses) && (length(mu_b1) != length(doses))) {
    stop("length(doses) must match length(mu_b1) and length(sigma_b1)")
  }
}
