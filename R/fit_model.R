fit_model <- function(
  model,
  data = NULL,
  doses,
  n_adapt,
  n_burn,
  n_iter,
  n_chains,
  jags_rng,
  jags_seed,
  silent
) {
  jset <- rjags_setup(
    silent = silent,
    jags_rng = model$jags_rng,
    jags_seed = model$jags_seed,
    n_chains = n_chains
  )
  jags_data <- get_jags_data(model, data, doses)
  jags_model <- get_jags_model(model, data)
  variable_names <- get_vnames(model)
  post_samples <- run_jags_model(
    jags_model = jags_model,
    jags_data = jags_data,
    jset = jset,
    variable_names = variable_names,
    n_adapt = n_adapt,
    n_burn = n_burn,
    n_iter = n_iter,
    n_chains = n_chains,
    silent = silent
  ) %>%
    add_mcmc_attributes(data, model, doses) %>%
    add_mcmc_class(model)
  return(post_samples)
}

get_progress_bar <- function(silent) {
  if (silent) {
    progress_bar <- "none"
  } else {
    progress_bar <- "text"
  }
  return(progress_bar)
}

get_default_jags_rng <- function(jags_rng) {
  if (is.null(jags_rng)) {
    jags_rng <- "base::Mersenne-Twister"
  }
  return(jags_rng)
}

expand_jags_rng <- function(jags_rng, n_chains) {
  if (length(jags_rng) == 1) {
    jags_rng <- rep(jags_rng, n_chains)
  } else if (length(jags_rng) != n_chains) {
    stop("jags_rng must have length 1 or n_chains.", call. = FALSE)
  }
  return(jags_rng)
}

check_seed_len <- function(jags_seed, n_chains) {
  if (!is.null(jags_seed)) {
    if (length(jags_seed) != n_chains) {
      stop(
        "jags_seed must be NULL or specify a seed for each chain",
        call. = FALSE
      )
    }
  }
}

get_jags_seed <- function(jags_seed, jags_rng, n_chains) {
  jags_inits <- NULL
  jags_inits <- list()
  check_seed_len(jags_seed, n_chains)
  for (i in 1:n_chains) {
    jags_inits[[i]] <- set_jags_seed(jags_seed[i], jags_rng[i])
  }
  return(jags_inits)
}

set_jags_seed <- function(jags_seed, jags_rng) {
  if (is.null(jags_seed)) {
    jags_seed <- sample(.Machine$integer.max, 1)
  }
  list(
    ".RNG.name" = jags_rng,
    ".RNG.seed" = jags_seed
  )
}

rjags_setup <- function(silent, jags_rng, jags_seed, n_chains) {
  progress_bar <- get_progress_bar(silent)
  jags_rng <- get_default_jags_rng(jags_rng) %>%
    expand_jags_rng(n_chains)
  jags_inits <- get_jags_seed(jags_seed, jags_rng, n_chains)
  if (!("glm" %in% rjags::list.modules()))
    suppressMessages(rjags::load.module("glm"))
  return(list(
    progress_bar = progress_bar,
    jags_inits = jags_inits
  ))
}

run_jags_model <- function(
  jags_model,
  jags_data,
  jset,
  variable_names,
  n_adapt,
  n_burn,
  n_iter,
  n_chains,
  silent
) {
  jmod <- rjags::jags.model(
    file = textConnection(jags_model),
    data = jags_data,
    inits = jset$jags_inits,
    n.adapt = n_adapt,
    n.chains = n_chains,
    quiet = silent
  )
  if (n_burn > 0) {
    stats::update(jmod, n.iter = n_burn, progress.bar = jset$progress_bar)
  }
  samps <- rjags::coda.samples(
    jmod,
    variable.names = variable_names,
    n.iter = n_iter,
    progress.bar = jset$progress_bar
  )
  return(samps)
}

add_mcmc_attributes <- function(samps, data, model, doses) {
  UseMethod("add_mcmc_attributes", model)
}

add_mcmc_attributes.default <- function(samps, data, model, doses) {
  common_attributes(samps, data, model, doses)
}

add_mcmc_attributes.dreamer_binary <- function(samps, data, model, doses) { #nolint
  samps <- common_attributes(samps, data, model, doses) %>%
    add_link_attr(model)
  return(samps)
}

add_mcmc_attributes.dreamer_beta <- function(samps, data, model, doses) { #nolint
  samps <- samps %>%
    common_attributes(data, model, doses)
  attr(samps, "scale") <- get_scale(model, data)
  return(samps)
}

add_mcmc_attributes.dreamer_beta_binary <- function(samps, data, model, doses) { #nolint
  samps <- samps %>%
    common_attributes(data, model, doses) %>%
    add_link_attr(model)
  attr(samps, "scale") <- get_scale(model, data)
  return(samps)
}

common_attributes <- function(samps, data, model, doses) {
  times <- NULL
  long_mod_attribute <- NULL
  if (!is.null(model$longitudinal)) {
    times <- sort(unique(data$time))
    long_mod_attribute <- class(model$longitudinal)
  }
  attr(samps, "doses") <- doses
  attr(samps, "times") <- times
  attr(samps, "longitudinal_model") <- long_mod_attribute
  attr(samps, "t_max") <- model$longitudinal$t_max
  return(samps)
}

add_link_attr <- function(x, model) {
  attr(x, "link") <- model$link
  return(x)
}

get_scale <- function(model, data) {
  if (is.null(model$scale) & !is.null(data)) {
    model$scale <- 1.2 * max(data$dose)
  } else if (is.null(model$scale) & is.null(data)) {
    stop("please specify scale (1.2 * (maximum investigational dose)?)")
  }
  return(model$scale)
}

add_mcmc_class <- function(post_samples, model) {
  UseMethod("add_mcmc_class", model)
}

add_mcmc_class.dreamer_linear <- function(post_samples, model) {
  class(post_samples) <- c(
    "dreamer_mcmc_linear",
    "dreamer_mcmc_continuous",
    "dreamer_mcmc",
    "dreamer",
    class(post_samples)
  )
  return(post_samples)
}

add_mcmc_class.dreamer_quad <- function(post_samples, model) {
  class(post_samples) <- c(
    "dreamer_mcmc_quad",
    "dreamer_mcmc_continuous",
    "dreamer_mcmc",
    "dreamer",
    class(post_samples)
  )
  return(post_samples)
}

add_mcmc_class.dreamer_loglinear <- function(post_samples, model) { #nolint
  class(post_samples) <- c(
    "dreamer_mcmc_loglinear",
    "dreamer_mcmc_continuous",
    "dreamer_mcmc",
    "dreamer",
    class(post_samples)
  )
  return(post_samples)
}

add_mcmc_class.dreamer_logquad <- function(post_samples, model) {
  class(post_samples) <- c(
    "dreamer_mcmc_logquad",
    "dreamer_mcmc_continuous",
    "dreamer_mcmc",
    "dreamer",
    class(post_samples)
  )
  return(post_samples)
}

add_mcmc_class.dreamer_emax <- function(post_samples, model) {
  class(post_samples) <- c(
    "dreamer_mcmc_emax",
    "dreamer_mcmc_continuous",
    "dreamer_mcmc",
    "dreamer",
    class(post_samples)
  )
  return(post_samples)
}

add_mcmc_class.dreamer_exp <- function(post_samples, model) {
  class(post_samples) <- c(
    "dreamer_mcmc_exp",
    "dreamer_mcmc_continuous",
    "dreamer_mcmc",
    "dreamer",
    class(post_samples)
  )
  return(post_samples)
}

add_mcmc_class.dreamer_beta <- function(post_samples, model) {
  class(post_samples) <- c(
    "dreamer_mcmc_beta",
    "dreamer_mcmc_continuous",
    "dreamer_mcmc",
    "dreamer",
    class(post_samples)
  )
  return(post_samples)
}

add_mcmc_class.dreamer_independent_continuous <- function(post_samples, model) { #nolint
  class(post_samples) <- c(
    "dreamer_mcmc_independent",
    "dreamer_mcmc_continuous",
    "dreamer_mcmc",
    "dreamer",
    class(post_samples)
  )
  return(post_samples)
}

add_mcmc_class.dreamer_linear_binary <- function(post_samples, model) { #nolint
  class(post_samples) <- c(
    "dreamer_mcmc_linear_binary",
    "dreamer_mcmc_binary",
    "dreamer_mcmc",
    "dreamer",
    class(post_samples)
  )
  return(post_samples)
}

add_mcmc_class.dreamer_quad_binary <- function(post_samples, model) { #nolint
  class(post_samples) <- c(
    "dreamer_mcmc_quad_binary",
    "dreamer_mcmc_binary",
    "dreamer_mcmc",
    "dreamer",
    class(post_samples)
  )
  return(post_samples)
}

add_mcmc_class.dreamer_loglinear_binary <- function(post_samples, model) { #nolint
  class(post_samples) <- c(
    "dreamer_mcmc_loglinear_binary",
    "dreamer_mcmc_binary",
    "dreamer_mcmc",
    "dreamer",
    class(post_samples)
  )
  return(post_samples)
}

add_mcmc_class.dreamer_logquad_binary <- function(post_samples, model) { #nolint
  class(post_samples) <- c(
    "dreamer_mcmc_logquad_binary",
    "dreamer_mcmc_binary",
    "dreamer_mcmc",
    "dreamer",
    class(post_samples)
  )
  return(post_samples)
}

add_mcmc_class.dreamer_emax_binary <- function(post_samples, model) { #nolint
  class(post_samples) <- c(
    "dreamer_mcmc_emax_binary",
    "dreamer_mcmc_binary",
    "dreamer_mcmc",
    "dreamer",
    class(post_samples)
  )
  return(post_samples)
}

add_mcmc_class.dreamer_exp_binary <- function(post_samples, model) { #nolint
  class(post_samples) <- c(
    "dreamer_mcmc_exp_binary",
    "dreamer_mcmc_binary",
    "dreamer_mcmc",
    "dreamer",
    class(post_samples)
  )
  return(post_samples)
}

add_mcmc_class.dreamer_beta_binary <- function(post_samples, model) { #nolint
  class(post_samples) <- c(
    "dreamer_mcmc_beta_binary",
    "dreamer_mcmc_binary",
    "dreamer_mcmc",
    "dreamer",
    class(post_samples)
  )
  return(post_samples)
}

add_mcmc_class.dreamer_independent_binary <- function(post_samples, model) { #nolint
  class(post_samples) <- c(
    "dreamer_mcmc_independent_binary",
    "dreamer_mcmc_binary",
    "dreamer_mcmc",
    "dreamer",
    class(post_samples)
  )
  return(post_samples)
}
