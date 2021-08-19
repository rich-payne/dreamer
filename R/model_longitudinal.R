#' @title Model Creation: Longitudinal Models
#' @name model_longitudinal
#' @description Assign hyperparameters and other values for longitudinal
#'   modeling.  The output of this function is intended to be used as
#'   the input to the `longitudinal` argument of the dose response model
#'   functions, e.g., `model_linear`.
#' @param t_max a scalar, typically indicating the latest observed time
#'   for subjects.  This will influence the interpretation of the
#'   parameters of each model.
#' @param mu_a,sigma_a,a_c1,b_c1,a_c2,b_c2 hyperparameters of the specified
#'   longitudinal model.  See below for parameterization.
#' @return A named list of the arguments in the function call.  The list has
#'   S3 classes assigned which are used internally within `dreamer_mcmc()`.
NULL

#' @rdname model_longitudinal
#' @section Longitudinal Linear:
#'   Let \eqn{f(d)} be a dose response model.  The expected value of the
#'   response, y, is:
#'   \deqn{E(y) = g(d, t)}
#'   \deqn{g(d, t) = a + (t / t_max) * f(d)}
#'   \deqn{a \sim N(mu_a, sigma_a)}
#' @export
model_longitudinal_linear <- function(
  mu_a,
  sigma_a,
  t_max
) {
  out <- list(
    t_max = t_max,
    mu_a = mu_a,
    sigma_a = sigma_a
  )
  class(out) <- c("dreamer_longitudinal_linear", "dreamer_longitudinal")
  return(out)
}

#' @rdname model_longitudinal
#' @section Longitudinal ITP:
#'   Let \eqn{f(d)} be a dose response model.  The expected value of the
#'   response, y, is:
#'   \deqn{E(y) = g(d, t)}
#'   \deqn{g(d, t) = a + f(d) * ((1 - exp(- c1 * t))/(1 - exp(- c1 * t_max)))}
#'   \deqn{a \sim N(mu_a, sigma_a)}
#'   \deqn{c1 \sim Uniform(a_c1, b_c1)}
#' @export
model_longitudinal_itp <- function(
  mu_a,
  sigma_a,
  a_c1 = 0,
  b_c1 = 1,
  t_max
) {
  out <- list(
    mu_a = mu_a,
    sigma_a = sigma_a,
    a_c1 = a_c1,
    b_c1 = b_c1,
    t_max = t_max
  )
  class(out) <- c("dreamer_longitudinal_itp", "dreamer_longitudinal")
  return(out)
}

#' @rdname model_longitudinal
#' @section Longitudinal IDP:
#'   Increasing-Decreasing-Plateau (IDP).
#'
#'   Let \eqn{f(d)} be a dose response model.  The expected value of the
#'   response, y, is:
#'   \deqn{E(y) = g(d, t)}
#'   \deqn{g(d, t) = a + f(d) * (((1 - exp(- c1 * t))/(1 - exp(- c1 * d1))) *
#'     I(t < d1) + (1 - gam * ((1 - exp(- c2 * (t - d1))) /
#'     (1 - exp(- c2 * (d2 - d1))))) *
#'     I(d1 <= t <= d2) + (1 - gam) * I(t > d2))}
#'   \deqn{a \sim N(mu_a, sigma_a)}
#'   \deqn{c1 \sim Uniform(a_c1, b_c1)}
#'   \deqn{c2 \sim Uniform(a_c2, b_c2)}
#'   \deqn{d1 \sim Uniform(0, t_max)}
#'   \deqn{d2 \sim Uniform(d1, t_max)}
#'   \deqn{gam \sim Uniform(0, 1)}
#' @export
model_longitudinal_idp <- function(
  mu_a,
  sigma_a,
  a_c1 = 0,
  b_c1 = 1,
  a_c2 = -1,
  b_c2 = 0,
  t_max
) {
  out <- list(
    mu_a = mu_a,
    sigma_a = sigma_a,
    a_c1 = a_c1,
    b_c1 = b_c1,
    a_c2 = a_c2,
    b_c2 = b_c2,
    t_max = t_max
  )
  class(out) <- c("dreamer_longitudinal_idp", "dreamer_longitudinal")
  return(out)
}
