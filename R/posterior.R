#' @title Posterior Quantities from Bayesian Model Averaging
#' @description Calculate posterior mean (and quantiles for specific
#'   doses for each MCMC iteration of the model.
#' @param x output from a call to \code{\link{dreamer_mcmc}}.
#' @param doses doses at which to estimate posterior quantities.
#' @param times a vector of times at which to calculate the posterior response
#'   (for longitudinal models only).
#' @param probs quantiles of the posterior to be calculated.
#' @param reference_dose the dose at which to adjust the posterior plot.
#'   Specifying
#'   a dose returns the plot of pr(trt_dose - trt_{reference_dose} | data).
#' @param return_samples logical indicating if the weighted raw
#'   MCMC samples from the Bayesian model averaging used to calculate the
#'   mean and quantiles should be returned.
#' @param predictive An integer.  If greater than 0, the return values will
#'   be from the predictive distribution of the mean of `predictive`
#'   observations.
#'   If 0 (default), the posterior on the dose response mean is returned.
#' @param iter an index on which iterations of the MCMC should be used
#'   in the calculations.  By default, all MCMC iterations are used.
#' @param return_stats logical indicating whether or not the posterior
#'   statistics should be calculated.
#' @return A named list with the following elements:
#'   * stats: a tibble the dose, posterior mean, and posterior quantiles.
#'   * samps: the weighted posterior samples.  Only returned if
#'     `return_samples = TRUE`.
#' @example man/examples/ex-posterior.R
#' @export
posterior <- function(
  x,
  doses,
  times,
  probs,
  reference_dose,
  predictive,
  return_samples,
  iter,
  return_stats
) {
  UseMethod("posterior", x)
}

#' @export
posterior.default <- function(
  x,
  doses,
  times,
  probs,
  reference_dose,
  predictive,
  return_samples,
  iter,
  return_stats
) {
  stop(
    "Class ",
    paste0(class(x), collapse = ", "),
    "not supported"
  )
}

get_post_mean_samps <- function(x, doses, times, iter = NULL) {
  if (is.null(iter)) {
    iter <- 1:(length(x) * nrow(x[[1]]))
  }
  if (is.null(times)) {
    the_grid <- data.frame(dose = doses)
    the_grid$time <- list(NULL)
    out <- expand.grid(iter = iter, dose = doses, KEEP.OUT.ATTRS = FALSE)
    out$time <- list(NULL)
  } else {
    the_grid <- expand.grid(dose = doses, time = times, KEEP.OUT.ATTRS = FALSE)
    out <- expand.grid(
      iter = iter,
      dose = doses,
      time = times,
      KEEP.OUT.ATTRS = FALSE
    )
  }
  y <- subset_mcmc(x, iter)
  vals <- apply(
    the_grid,
    1,
    function(xx) dreamer_mean(x = x, y = y, xx[[1]], xx[[2]], iter)
  )
  out <- out %>%
    dplyr::mutate(mean_response = c(vals))
  return(out)
}

#' @describeIn posterior posterior summary for linear model.
#' @export
posterior.dreamer_mcmc <- function(
  x,
  doses = attr(x, "doses"),
  times = attr(x, "times"),
  probs = c(.025, .975),
  reference_dose = NULL,
  predictive = 0,
  return_samples = FALSE,
  iter = NULL,
  return_stats = TRUE
) {
  assert_reference_dose(reference_dose)
  mcmc_samps <- subset_mcmc(x, iter)
  original_doses <- doses
  doses <- sort(unique(c(reference_dose, doses)))
  samps <- get_post_mean_samps(x = x, doses = doses, times = times, iter = iter)
  output <- summarize_samples(
    x = x,
    mcmc_samps = mcmc_samps,
    samps = samps,
    probs = probs,
    doses = doses,
    original_doses = original_doses,
    return_samples = return_samples,
    reference_dose = reference_dose,
    return_stats = return_stats,
    predictive = predictive
  )
  return(output)
}

#' @describeIn posterior posterior summary for independent model.
#' @export
posterior.dreamer_mcmc_independent <- function( #nolint
  x,
  doses = attr(x, "doses"),
  times = attr(x, "times"),
  probs = c(.025, .975),
  reference_dose = NULL,
  predictive = 0,
  return_samples = FALSE,
  iter = NULL,
  return_stats = TRUE
) {
  all_doses <- attributes(x)$doses
  original_doses <- doses
  doses <- sort(unique(c(reference_dose, doses)))
  check_independent_doses(doses = doses, all_doses = all_doses)
  mcmc_samps <- as.matrix(x)
  if (!is.null(iter)) {
    mcmc_samps <- mcmc_samps[iter, , drop = FALSE]
  }
  samps <- get_post_mean_samps(x = x, doses = doses, times = times, iter = iter)
  output <- summarize_samples(
    x = x,
    mcmc_samps = mcmc_samps,
    samps = samps,
    probs = probs,
    doses = doses,
    original_doses = original_doses,
    return_samples = return_samples,
    reference_dose = reference_dose,
    return_stats = return_stats,
    predictive = predictive,
    response_type = attr(x, "response_type")
  )
  return(output)
}

check_independent_doses <- function(doses, all_doses) {
  if (!is.null(doses)) {
    if (!isTRUE(all(doses %in% all_doses))) {
      stop(
        "All doses must be in the original data for an independent model.",
        call. = FALSE
      )
    }
  }
}

#' @describeIn posterior posterior summary for Bayesian model averaging fit.
#' @export
posterior.dreamer_bma <- function(
  x,
  doses = x$doses,
  times = x$times,
  probs = c(.025, .975),
  reference_dose = NULL,
  predictive = 0,
  return_samples = FALSE,
  iter = NULL,
  return_stats = TRUE
) {
  mcmc_index <- get_mcmc_index(x)
  if (is.null(iter)) {
    iter <- seq_len(length(attr(x, "model_index")))
  }
  model_index <- attr(x, "model_index")[iter]
  u_model_index <- unique(model_index)
  post_samps <- purrr::map_df(
    u_model_index,
    function(u_model_index) {
      ind <- which(model_index == u_model_index)
      ind_iter <- iter[ind]
      out <- posterior(
        x = x[[mcmc_index[u_model_index]]],
        doses = doses,
        times = times,
        probs = probs,
        reference_dose = reference_dose,
        predictive = predictive,
        return_samples = TRUE,
        iter = ind_iter,
        return_stats = FALSE
      )$samps
    }
  )
  output <- summarize_samples(
    x = NULL, # not needed here
    samps = post_samps,
    probs = probs,
    doses = doses,
    original_doses = x$doses,
    return_samples = return_samples,
    reference_dose = NULL, # already adjusted above
    predictive = 0, # already adjusted above
    return_stats = return_stats,
    response_type = attr(x, "response_type")
  )
  return(output)
}

get_mcmc_index <- function(x) {
  vapply(
    x,
    function(y) any(grepl("mcmc_", class(y))),
    logical(1)
  ) %>%
    which()
}

logit <- function(x) {
  log(x / (1 - x))
}

probit <- function(x) {
  qnorm(x)
}

ilogit <- function(x) {
  1 / (1 + exp(-x))
}

iprobit <- function(x) {
  pnorm(x)
}

summarize_samples <- function(
  x,
  mcmc_samps,
  samps,
  probs,
  doses,
  original_doses,
  return_samples,
  reference_dose,
  return_stats = TRUE,
  response_type,
  predictive
) {
  samps <- dplyr::select(
    samps,
    tidyselect::vars_select_helpers$where((~ !is.list(.x)))
  )
  if (inherits(x, "dreamer_mcmc_continuous")) {
    samps <- summarize_samples_continuous(
      mcmc_samps = mcmc_samps,
      samps = samps,
      reference_dose = reference_dose,
      original_doses = original_doses,
      predictive = predictive
    )
  } else if (inherits(x, "dreamer_mcmc_binary")) {
    samps <- summarize_samples_binary(
      mcmc_samps = mcmc_samps,
      samps = samps,
      reference_dose = reference_dose,
      original_doses = original_doses,
      predictive = predictive
    )
  }
  samps <- samps %>%
    add_model_name(x) %>%
    add_reference_dose(reference_dose)
  stats <- NULL
  if (return_stats) {
    stats <- samps %>%
      dplyr::group_by(across(any_of(c("dose", "reference_dose", "time")))) %>%
      dplyr::summarize(
        mean = mean(.data$mean_response),
        qtiles = list(quantile_custom(.data$mean_response, !!probs))
      ) %>%
      dplyr::ungroup() %>%
      tidyr::unnest_wider(col = .data$qtiles)
  }
  output <- list(stats = stats)
  if (return_samples) {
    output$samps <- samps
  }
  return(output)
}

add_model_name <- function(samps, x) {
  if (!is.null(x)) {
    samps <- samps %>%
      dplyr::mutate(model = attr(x, "model_name"))
  }
  samps
}

add_reference_dose <- function(x, reference_dose) {
  if (!is.null(reference_dose)) {
    x <- x %>%
      dplyr::mutate(reference_dose = reference_dose)
  }
  x
}

quantile_custom <- function(x, probs, ...) {
  q <- quantile(x, probs = probs, names = FALSE, ...)
  names(q) <- vapply(
    100 * probs,
    sprintf,
    fmt = "%.2f%%",
    character(1)
  )
  q
}

make_predictive_continuous <- function(
  mcmc_samps,
  samps,
  predictive,
  reference_dose
) {
  if (!is.null(reference_dose)) {
    predictive <- predictive / 2
  }
  if (predictive > 0) {
    samps <- samps %>%
      dplyr::group_by(.data$iter) %>%
      tidyr::nest() %>%
      dplyr::ungroup() %>%
      dplyr::mutate(sigma = mcmc_samps[, "sigma"]) %>%
      tidyr::unnest("data") %>%
      dplyr::mutate(
        mean_response = .data$mean_response +
          rnorm(n(), 0, .data$sigma / sqrt(!!predictive))
      )
  }
  return(samps)
}

make_predictive_binary <- function(samps, predictive) {
  if (predictive > 0) {
    samps <- samps %>%
      dplyr::mutate(
        mean_response = rbinom(n(), !!predictive, .data$mean_response) /
          !!predictive
      )
  }
  return(samps)
}

make_reference_dose <- function(samps, reference_dose, original_doses) {
  add_adjustment(samps, reference_dose) %>%
    dplyr::mutate(
      mean_response = .data$mean_response - .data$mean_response_adj
    ) %>%
    dplyr::select(-.data$reference_dose, -.data$mean_response_adj) %>%
    dplyr::filter(.data$dose %in% !!original_doses)
}

add_adjustment <- function(samps, reference_dose) {
  by_var <- "iter"
  if (rlang::has_name(samps, "time")) {
    by_var <- c(by_var, "time")
  }
  samps_adj <- samps %>%
    dplyr::filter(.data$dose == !!reference_dose) %>%
    dplyr::select(reference_dose = .data$dose, everything())
  samps <- samps %>% dplyr::left_join(
    samps_adj,
    by = by_var,
    suffix = c("", "_adj")
  )
  return(samps)
}

summarize_samples_continuous <- function(
  mcmc_samps,
  samps,
  reference_dose,
  original_doses,
  predictive
) {
  if (!is.null(reference_dose)) {
    samps <- make_reference_dose(samps, reference_dose, original_doses)
  }
  if (predictive > 0) {
    samps <- make_predictive_continuous(
      mcmc_samps = mcmc_samps,
      samps = samps,
      predictive = predictive,
      reference_dose = reference_dose
    )
  }
  samps
}

summarize_samples_binary <- function(
  mcmc_samps,
  samps,
  reference_dose,
  original_doses,
  predictive
) {
  if (!is.null(reference_dose) & predictive == 0) {
    samps <- make_reference_dose(samps, reference_dose, original_doses)
  } else if (!is.null(reference_dose) & predictive > 0) {
    samps <- add_adjustment(samps, reference_dose) %>%
      dplyr::mutate(
        mean_response = (
            rbinom(n(), !!predictive, .data$mean_response) -
            rbinom(n(), !!predictive, .data$mean_response_adj)
          ) / !!predictive
      ) %>%
      dplyr::filter(.data$dose %in% !!original_doses) %>%
      dplyr::select(-.data$reference_dose, -.data$mean_response_adj)
  } else if (is.null(reference_dose)) {
    samps <- make_predictive_binary(samps, predictive)
  }
  return(samps)
}
