#' @title Bayesian Model Averaging of Dose Response Models
#' @description This function performs Bayesian model averaging with a selection
#'   of dose response models.  See \link[dreamer]{model} for all possible
#'   models.
#' @param data a dataframe with column names of "dose" and "response" for
#'   individual patient data.  Optional columns "n" and "sample_var" can
#'   be specified if aggregate data is supplied, but it is recommended that
#'   patient-level data be supplied where possible for continuous models, as the
#'   posterior weights differ if aggregated data is used.
#'   For aggregated continuous
#'   data, "response" should be the average of "n" subjects with a sample
#'   variance of "sample_var".  For aggregated binary data, "response" should
#'   be the number of successes, "n" should be the total number of subjects
#'   (the "sample_var" column is irrelevant in binary cases and is ignored).
#' @param ... model definitions created using the model creation functions
#'   in \link[dreamer]{model}.  If arguments are named, the names are
#'   retained in the return values.
#' @param n_adapt the number of MCMC iterations to tune the MCMC algorithm.
#' @param n_burn the number of burn-in MCMC samples.
#' @param n_iter the number of MCMC samples to collect after tuning and burn-in.
#' @param n_chains the number of separate, independent, MCMC chains to run.
#' @param silent logical indicating if MCMC progress bars should be suppressed.
#' @param convergence_warn logical (default `TRUE`) indicating if the
#'   Gelman-Rubin diagnostics should be run to detect convergence issues.
#'   Warnings are thrown if the upper bound of the Gelman-Rubin statistic
#'   is greater than 1.1.
#' @return A named list with S3 class "dreamer_bma" and "dreamer".  The list
#'   contains the following fields:
#'   * doses: a vector of the unique ordered doses in the data.
#'   * times: a vector of the unique ordered times in the data.
#'   * w_prior: a named vector with the prior probabilities of each model.
#'   * w_post: a named vector with the posterior probabilities of each model.
#'   * The individual MCMC fits for each model.
#' @details The Bayesian model averaging approach uses data, multiple models,
#'   priors on each model's parameters, and a prior weight for each model.
#'   Using these inputs, each model is fit independently, and the output from
#'   the models is used to calculate posterior weights for each model.
#'   See Gould (2018) for details.
#' @section References:
#'   Gould, A. L. (2019).  BMA-Mod: A Bayesian model averaging strategy for
#'     determining dose-response relationships in the presence of model
#'     uncertainty. Biometrical Journal, 61(5), 1141-1159.
#' @example man/examples/ex-dreamer_mcmc.R
#' @export
dreamer_mcmc <- function( #nolint
  data,
  ...,
  n_adapt = 1e3,
  n_burn = 1e3,
  n_iter = 1e4,
  n_chains = 4,
  silent = FALSE,
  convergence_warn = TRUE
) {
  jags_modules <- rjags::list.modules()
  on.exit(restore_jags_modules(jags_modules))
  load_jags_modules()
  mods <- list(...)
  assert_dreamer_dots(mods)
  assert_independent_dots(mods)
  all_dots_binary <- assert_binary_dots(mods)
  check_data(data, all_dots_binary)
  is_long <- check_longitudinal(mods, data)
  if (!is.null(data)) {
    if (!all_dots_binary) {
      check_aggregate(data)
    }
  }
  doses <- get_doses(data)
  all_models <- name_models(mods)
  mcmc_list <- list()
  for (i in seq_along(all_models)) {
    start_time <- Sys.time()
    mcmc_start_msg(names(all_models)[i], start_time, silent)
    mcmc_list[[i]] <- fit_model(
      all_models[[i]],
      data = data,
      doses = doses,
      n_adapt = n_adapt,
      n_burn = n_burn,
      n_iter = n_iter,
      n_chains = n_chains,
      silent = silent,
      binary = all_dots_binary
    )
    end_time <- Sys.time()
    mcmc_end_msg(start_time, end_time, silent)
  }
  names(mcmc_list) <- names(all_models)

  w_prior <- get_w_prior(all_models)
  assert_w_prior(w_prior)
  weight_data <- prep_weight_data(data, all_dots_binary, is_long)
  w <- post_weights(weight_data, mcmc_list, w_prior)

  model_index <- get_model_index(mcmc_list, w, n_iter, n_chains)
  times <- if (is_long) sort(unique(data$time)) else NULL
  model_names <- names(all_models)

  final_output <- c(
    list(
      doses = doses,
      times = times,
      w_prior = w$w_prior,
      w_post = w$w_post
    ),
    mcmc_list
  ) %>%
    add_attributes(
      all_dots_binary,
      is_long,
      doses,
      times,
      model_index,
      model_names
    )

  class(final_output) <- c("dreamer_bma", "dreamer_mcmc")
  if (convergence_warn) convergence_warnings(x = final_output)
  return(final_output)
}

load_jags_modules <- function() {
  rjags::load.module("glm", quiet = TRUE)
}

restore_jags_modules <- function(original_modules) {
  current_modules <- rjags::list.modules()
  ind <- which(!(current_modules %in% original_modules))
  if (length(ind) > 0)
    purrr::walk(current_modules[ind], ~ rjags::unload.module(.x, quiet = TRUE))
}

get_doses <- function(data) {
  sort(unique(data$dose))
}

add_attributes <- function(
  final_output,
  all_dots_binary,
  is_long,
  doses,
  times,
  model_index,
  model_names
) {
  attr(final_output, "response_type") <- get_response_type(all_dots_binary)
  if (is_long) {
    attr(final_output, "longitudinal_model") <- TRUE
  } else {
    attr(final_output, "longitudinal_model") <- NULL
  }
  attr(final_output, "doses") <- doses
  attr(final_output, "times") <- times
  attr(final_output, "model_index") <- model_index
  attr(final_output, "model_names") <- model_names
  mcmc_index <- get_mcmc_index(final_output)
  final_output[mcmc_index] <- purrr::map2(
    final_output[mcmc_index],
    model_names,
    add_model_name_attr
  )
  return(final_output)
}

add_model_name_attr <- function(x, model_name) {
  attr(x, "model_name") <- model_name
  x
}

get_model_index <- function(mcmc_list, w, n_iter, n_chains) {
  names_all_post <- names(mcmc_list)
  names_index <- vapply(
    names(w$w_post),
    function(y) which(y == names_all_post),
    integer(1)
  )
  n_mcmc <- n_iter * n_chains
  model_index <- base::sample(
    x = names_index,
    size = n_mcmc,
    replace = TRUE,
    prob = w$w_post
  )
  return(model_index)
}

post_weights <- function(data, mcmc_list, w_prior) {
  if (!is.null(data)) {
    w <- dreamer_post_weights(mcmc_list, w_prior, data)
  } else {
    w <- list(w_prior = w_prior, w_post = w_prior)
  }
  return(w)
}

mcmc_start_msg <- function(model_name, start_time, silent) {
  if (!silent) {
    message(model_name)
    message(paste0("start : ", format(start_time, "%Y-%m-%d %H:%M:%OS3")))
  }
  return(start_time)
}

mcmc_end_msg <- function(start_time, end_time, silent) {
  if (!silent) {
    time_diff <- end_time - start_time
    time_unit <- attributes(time_diff)$units
    class(time_diff) <- NULL
    message(paste0("finish: ", format(end_time, "%Y-%m-%d %H:%M:%OS3")))
    message(
      paste0(
        "total : ", format(round(time_diff, 2), nsmall = 2),
        " ", time_unit, "\n"
      )
    )
  }
}

assert_w_prior <- function(w_prior) {
  w_total <- sum(w_prior)
  if (!isTRUE(all.equal(w_total, 1))) {
    stop(
      "Sum of w_prior for all models must add to 1: ", w_total,
      call. = FALSE
    )
  }
}

get_w_prior <- function(all_models) {
  vapply(
    all_models,
    function(model) model$w_prior,
    FUN.VALUE = numeric(1)
  )
}

assert_dreamer_dots <- function(mods) {
  all_dots_dreamer <- vapply(
    mods,
    function(model) {
      any(inherits(model, c("dreamer_continuous", "dreamer_binary")))
    },
    logical(1)
  ) %>%
    all()
  if (!all_dots_dreamer) {
    stop("All '...' arguments must be dreamer models")
  }
  return(all_dots_dreamer)
}

assert_independent_dots <- function(mods) {
  any_dots_independent <- vapply(
    mods,
    function(model) inherits(model, "dreamer_independent"),
    logical(1)
  ) %>%
    any()
  if (length(mods) > 1 && any_dots_independent) {
    stop("Independent model cannot be used in model averaging", call. = FALSE)
  }
}

assert_binary_dots <- function(mods) {
  all_dots_binary <- vapply(
    mods,
    function(model) inherits(model, "dreamer_binary"),
    logical(1)
  ) %>%
    all()
  any_dots_binary <- vapply(
    mods,
    function(model) inherits(model, "dreamer_binary"),
    logical(1)
  ) %>%
    any()
  if (any_dots_binary && !all_dots_binary) {
    stop("All models must be the same type: continuous/binary.", call. = FALSE)
  }
  return(all_dots_binary)
}

check_data <- function(dat, binary) {
  if (is.null(dat)) {
    return(NULL)
  }
  if (!all(c("dose", "response") %in% colnames(dat))) {
    rlang::abort(
      "data must have columns \"dose\" and \"response\".",
      class = "dreamer"
    )
  }
  if (!is.numeric(dat$response)) {
    rlang::abort(
      "Column \"response\" in data must be numeric.",
      class = "dreamer"
    )
  }
  if (!is.numeric(dat$dose)) {
    rlang::abort(
      "Column \"dose\" in data must be numeric.",
      class = "dreamer"
    )
  }
  if (
    binary &&
    !all(dat$response %in% c(0L, 1L)) &&
    !rlang::has_name(dat, "n")
  ) {
    rlang::abort(
      "Column \"response\" must contain only zeros and ones unless \"n\" is specified.", # nolint
      class = "dreamer"
    )
  }
}

check_aggregate <- function(dat) {
  cnames <- colnames(dat)
  if (any(c("sample_var", "n") %in% cnames)) {
    rlang::inform(
      paste0(
        "Specifying 'sample_var' and 'n' will give different posterior weights",
        " than non-aggregated data.\n",
        "Consider using non-aggregated data."
      ),
      class = "dreamer"
    )
  }
}

convergence_warnings <- function(x, ...) {
  summaries <- summary(x)
  bad_gelman <- dplyr::filter(summaries$summary, .data$gelman_upper > 1.1) %>%
    dplyr::select(
      .data$model, .data$param, .data$gelman_point, .data$gelman_upper
    )
  throw_convergence_warn(bad_gelman)
}

throw_convergence_warn <- function(bad_gelman) {
  if (nrow(bad_gelman) > 0) {
    begin_wrn <- paste0(
      "Convergence issues detected.\n",
      "Upper CI of Gelman-Rubin above 1.1 in the following:\n\n"
    )
    specific_warnings <- paste0(
      bad_gelman$model,
      ", ",
      bad_gelman$param,
      ": ",
      format(bad_gelman$gelman_upper, digits = 3),
      collapse = "\n"
    )
    warning(c(begin_wrn, specific_warnings), call. = FALSE)
  }
}

name_models <- function(model_list) {
  for (i in seq_along(model_list)) {
    if (
      isTRUE(names(model_list)[i] == "") ||
      is.null(names(model_list)[i]) ||
      is.na(names(model_list)[i])
    ) {
      names(model_list)[i] <- paste0("model_", i, "_", class(model_list[[i]]))
    }
  }
  if (any(duplicated(names(model_list)))) {
    stop("Duplicate model names are not allowed.", call. = FALSE)
  }
  return(model_list)
}

prep_weight_data <- function(data, binary, is_longitudinal) {
  if (binary) {
    data <- aggregate_binary(data, is_longitudinal)
  }
  return(data)
}

dreamer_post_weights <- function(all_mcmc, w_prior, data) {
  etas <- vapply(
    all_mcmc,
    function(x, data) epll(x, data),
    data = data,
    FUN.VALUE = numeric(1)
  )
  w_prior <- w_prior[names(etas)]
  assert_names(etas, w_prior, "Naming issue with etas/w_prior.")
  # log-sum-exp trick for numeric stability
  a_s <- etas + log(w_prior)
  a <- max(a_s)
  w_post <- exp(a_s - (a + log(sum(exp(a_s - a)))))
  names(w_post) <- names(all_mcmc)
  return(list(w_prior = w_prior, w_post = w_post))
}

assert_names <- function(x, y, msg) {
  if (!all(names(x) == names(y))) {
    stop(msg, call. = FALSE)
  }
}

check_binary <- function(data) {
  all_non_neg <- all(data$response >= 0)
  all_0_or_1 <- all(data$response %in% c(0, 1))
  n_null <- !rlang::has_name(data, "n")
  if (!all_non_neg) {
    stop("All values of data$response must be non-negative.", call. = FALSE)
  }
  if (!all_0_or_1 && n_null) {
    stop(
      "All values of data$response must be 0 or 1",
      " unless n is specified in data.",
      call. = FALSE
    )
  }
}

check_longitudinal <- function(mods, data) {
  is_long <- vapply(
    mods,
    function(x) !is.null(x$longitudinal),
    logical(1)
  )
  if (any(is_long) && !all(is_long)) {
    stop("All models must be longitudinal or non-longitudinal", call. = FALSE)
  }
  if (all(is_long) && (!rlang::has_name(data, "time") && !is.null(data))) {
    stop("data must have column 'time' if longitudinal models are specified.")
  }
  return(all(is_long))
}
