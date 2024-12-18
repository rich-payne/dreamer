# expected predicted log-likelihodd (epll)
epll <- function(x, data) {
  binary <- inherits(x, "dreamer_mcmc_binary")
  sum_log_probs <- get_sum_log_probs(x = x, data = data, binary = binary)
  q <- get_n_params(x)
  return(get_epll(sum_log_probs, q, n = nrow(data)))
}

get_n_params_longitudinal <- function(x) {
  long_mod <- attr(x, "longitudinal_model")
  if (is.null(long_mod)) {
    q <- 0
  } else if ("dreamer_longitudinal_itp" %in% long_mod) {
    q <- 1 + 1
  } else if ("dreamer_longitudinal_idp" %in% long_mod) {
    q <- 5 + 1
  } else if ("dreamer_longitudinal_linear" %in% long_mod) {
    q <- 0 + 1
  } else {
    stop(
      paste0(
        "Longitudinal model with class(es) '",
        paste0(long_mod, collapse = ", "),
        "' not yet supported."
      ),
      call. = FALSE
    )
  }
  return(q)
}

get_n_params <- function(x) {
  UseMethod("get_n_params", x)
}

#' @export
get_n_params.dreamer_mcmc_linear <- function(x) { #nolint
  q <- 3 + get_n_params_longitudinal(x)
  return(q)
}

#' @export
get_n_params.dreamer_mcmc_quad <- function(x) {
  q <- 4 + get_n_params_longitudinal(x)
  return(q)
}

#' @export
get_n_params.dreamer_mcmc_loglinear <- function(x) { #nolint
  q <- 3 + get_n_params_longitudinal(x)
  return(q)
}

#' @export
get_n_params.dreamer_mcmc_logquad <- function(x) { #nolint
  q <- 4 + get_n_params_longitudinal(x)
  return(q)
}

#' @export
get_n_params.dreamer_mcmc_emax <- function(x) {
  q <- ifelse(all(as.matrix(x)[, "b4"] == 1), 4, 5) +
    get_n_params_longitudinal(x)
  return(q)
}

#' @export
get_n_params.dreamer_mcmc_exp <- function(x) {
  q <- 4 + get_n_params_longitudinal(x)
  return(q)
}

#' @export
get_n_params.dreamer_mcmc_beta <- function(x) {
  q <- 5 + get_n_params_longitudinal(x)
  return(q)
}

#' @export
get_n_params.dreamer_mcmc_independent <- function(x) { #nolint
  # this function is not relevant when independent models are
  #   restricted to be fit in isolation
  return(1)
}

# BINARY
#' @export
get_n_params.dreamer_mcmc_linear_binary <- function(x) { #nolint
  q <- 2 + get_n_params_longitudinal(x)
  return(q)
}

#' @export
get_n_params.dreamer_mcmc_quad_binary <- function(x) { #nolint
  q <- 3 + get_n_params_longitudinal(x)
  return(q)
}

#' @export
get_n_params.dreamer_mcmc_loglinear_binary <- function(x) { #nolint
  q <- 2 + get_n_params_longitudinal(x)
  return(q)
}

#' @export
get_n_params.dreamer_mcmc_logquad_binary <- function(x) { #nolint
  q <- 3 + get_n_params_longitudinal(x)
  return(q)
}

#' @export
get_n_params.dreamer_mcmc_emax_binary <- function(x) { #nolint
  q <- ifelse(all(as.matrix(x)[, "b4"] == 1), 3, 4) +
    get_n_params_longitudinal(x)
  return(q)
}

#' @export
get_n_params.dreamer_mcmc_exp_binary <- function(x) { #nolint
  q <- 3 + get_n_params_longitudinal(x)
  return(q)
}

#' @export
get_n_params.dreamer_mcmc_beta_binary <- function(x) { #nolint
  q <- 4 + get_n_params_longitudinal(x)
  return(q)
}

#' @export
get_n_params.dreamer_mcmc_independent_binary <- function(x) { #nolint
  # this function is not relevant when independent models are
  #   restricted to be fit in isolation
  return(1)
}
# q is number of parameters
get_epll <- function(sum_log_probs, q, n) {
  eta <- sum_log_probs - q / 2
  return(eta)
}

get_sum_log_probs <- function(x, data, binary = FALSE) {
  post <- posterior(x, return_samples = TRUE)
  post_samps <- post$samps
  if (!binary) {
    sigmas <- as.matrix(x)[, "sigma"]
    post_sigma <- data.frame(
      iter = seq_len(length(sigmas)),
      sigma = sigmas
    )
    post_samps <- dplyr::left_join(post_samps, post_sigma, by = "iter")
  }
  if (!is.null(attr(x, "longitudinal_model"))) {
    obs_ind <- purrr::map2_int(
      data$dose,
      data$time,
      function(x, y, post) {
        which(x == post$stats$dose & y == post$stats$time)
      },
      post = post
    )
    mcmc_means <- purrr::pmap(
      post$stats,
      function(post_samps, dose, time, ...) {
        dplyr::filter(post_samps, dose == !!dose, time == !!time)
      },
      post_samps = post_samps
    )
  } else {
    obs_ind <- purrr::map_int(
      data$dose,
      function(x, post) {
        which(x == post$stats$dose)
      },
      post = post
    )
    mcmc_means <- purrr::pmap(
      post$stats,
      function(post_samps, dose, ...) {
        dplyr::filter(post_samps, dose == !!dose)
      },
      post_samps = post_samps
    )
  }
  if (!binary) {
    sum_log_probs <- purrr::map2_dbl(
      obs_ind,
      data$response,
      function(x, y) {
        vals <- dnorm(
          y,
          mcmc_means[[x]][["mean_response"]],
          mcmc_means[[x]][["sigma"]]
        ) %>%
          mean()
        return(vals)
      }
    ) %>%
      log() %>%
      sum()
  } else {
    arg_list <- data.frame(
      obs_ind = obs_ind,
      response = data$response,
      n = data$n
    )
    sum_log_probs <- purrr::pmap_dbl(
      arg_list,
      function(obs_ind, response, n) {
        vals <- dbinom(
          x = response,
          size = n,
          prob = mcmc_means[[obs_ind]][["mean_response"]]
        ) %>%
          mean()
        return(vals)
      }
    ) %>%
      log() %>%
      sum()
  }
  return(sum_log_probs)
}
