base_dreamer_plot_error_bar <- function(
  x,
  doses,
  time,
  probs,
  predictive,
  width,
  reference_dose = NULL,
  p = NULL,
  color = "black",
  linetype = 1
) {
  check_plot(doses = doses, probs = probs)
  post <- posterior(
    x = x,
    doses = doses,
    times = time,
    probs = probs,
    predictive = predictive,
    return_samples = FALSE,
    reference_dose = reference_dose
  )
  dat <- data.frame(
    dose = post$stats$dose,
    post_mean = post$stats$mean,
    lb = post$stats[[sprintf("%.2f%%", 100 * probs[1])]],
    ub = post$stats[[sprintf("%.2f%%", 100 * probs[2])]]
  )
  if (is.null(p)) {
    p <- ggplot(dat, aes_string(x = "dose"))
  }
  p <- p +
    geom_errorbar(
    data = dat,
    mapping = aes_string(x = "dose", ymin = "lb", ymax = "ub", width = "width"),
    colour = color,
    linetype = linetype
  ) +
    scale_x_continuous(breaks = doses) +
    labs(x = "Dose", y = "Response") +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(title = get_title_error_bars(probs))
  return(p)
}

base_dreamer_plot <- function(
  x,
  doses,
  time,
  probs,
  reference_dose = NULL,
  p = NULL,
  color = "cornflowerblue"
) {
  post <- posterior(
    x = x,
    doses = doses,
    times = time,
    probs = probs,
    return_samples = FALSE,
    reference_dose = reference_dose
  )
  dat <- data.frame(
    dose = post$stats$dose,
    post_mean = post$stats$mean,
    lb = post$stats[[sprintf("%.2f%%", 100 * probs[1])]],
    ub = post$stats[[sprintf("%.2f%%", 100 * probs[2])]]
  )
  the_title <- get_title(probs, "Posterior") %>%
    add_time_to_title(time)
  p <- p +
    geom_line(
      data = dat,
      mapping = aes_string(x = "dose", y = "post_mean"),
      colour = color
    ) +
    geom_line(
      data = dat,
      mapping = aes_string(x = "dose", y = "lb"),
      colour = color,
      linetype = 2
    ) +
    geom_line(
      data = dat,
      mapping = aes_string(x = "dose", y = "ub"),
      colour = color,
      linetype = 2
    ) +
    labs(title = the_title)
  return(p)
}

add_time_to_title <- function(the_title, time) {
  if (length(time) == 1) {
    the_title <- paste0(the_title, ", time = ", time)
  }
  the_title
}

base_dreamer_predictive_plot <- function(
  x,
  doses,
  time,
  probs,
  predictive,
  reference_dose = NULL,
  p = NULL,
  color = "forestgreen"
) {
  post <- posterior(
    x = x,
    doses = doses,
    times = time,
    probs = probs,
    predictive = predictive,
    return_samples = FALSE,
    reference_dose = reference_dose
  )
  dat <- data.frame(
    dose = post$stats$dose,
    post_mean = post$stats$mean,
    lb = post$stats[[sprintf("%.2f%%", 100 * probs[1])]],
    ub = post$stats[[sprintf("%.2f%%", 100 * probs[2])]]
  )
  p <- p +
    geom_line(
      data = dat,
      mapping = aes_string(x = "dose", y = "lb"),
      colour = color,
      linetype = 3
    ) +
    geom_line(
      data = dat,
      mapping = aes_string(x = "dose", y = "ub"),
      colour = color,
      linetype = 3
    )
  return(p)
}

#' @name dreamerplot
#' @title Posterior Plot of Bayesian Model Averaging
#' @description Plots the posterior mean and quantiles over the dose range and
#'   plots error bars at the observed doses. If the `data` argument is
#'   specified, the observed means at each dose are also plotted.
#' @param x output from a call to \code{\link{dreamer_mcmc}}.
#' @param doses a vector of doses at which to plot the dose response curve.
#' @param probs quantiles of the posterior to be calculated.
#' @param times a vector of the times at which to plot the posterior (for
#'   longitudinal models only).
#' @inheritParams dreamer_mcmc
#' @param n_smooth the number of points to calculate the smooth dose response
#'   interpolation.  Must be sufficiently high to accurately depict the
#'   dose response curve.
#' @param predictive the size of sample for which to plot posterior predictive
#'   intervals for the mean.
#' @param width the width of the error bars.
#' @param reference_dose the dose at which to adjust the posterior plot.
#'   Specifying
#'   a dose returns the plot of pr(trt_dose - trt_{reference_dose} | data).
#' @return Returns the ggplot object.
#' @example man/examples/ex-plot.R
#' @export
plot.dreamer_mcmc <- function(
  x,
  doses = attr(x, "doses"),
  times = attr(x, "times"),
  probs = c(.025, .975),
  data = NULL,
  n_smooth = 50,
  predictive = 0,
  width = bar_width(doses),
  reference_dose = NULL,
  ...
) {
  check_no_dots("plot.dreamer()", ...)
  assert_data_reference_dose(data, reference_dose)
  times <- get_time(x, times, max_length = Inf)
  force(width)
  force(doses)
  if (!is.null(attr(x, "longitudinal_model")) & length(times) > 1) {
    p <- plot_longitudinal(
      x = x,
      doses = doses,
      times = times,
      probs = probs,
      n_smooth = n_smooth,
      predictive = predictive,
      reference_dose = reference_dose,
      data = data
    )
    return(p)
  }
  post <- posterior(
    x = x,
    doses = doses,
    times = times,
    probs = probs,
    return_samples = FALSE,
    reference_dose = reference_dose
  )
  dat <- dplyr::tibble(
    dose = post$stats$dose,
    post_mean = post$stats$mean,
    lb = post$stats[[sprintf("%.2f%%", 100 * probs[1])]],
    ub = post$stats[[sprintf("%.2f%%", 100 * probs[2])]]
  )
  p <- base_dreamer_plot_error_bar(
    x = x,
    doses = doses,
    time = times,
    probs = probs,
    reference_dose = reference_dose,
    predictive = 0, # add predictive in later if needed
    width = width
  )
  # continuous interpolation
  any_indep <- any_independent(x)
  range_of_doses <- seq(
    from = min(doses),
    to = max(doses),
    length.out = n_smooth
  )
  if (n_smooth > 0 & !any_indep) {
    p <- base_dreamer_plot(
      x = x,
      doses = range_of_doses,
      time = times,
      probs = probs,
      reference_dose = reference_dose,
      p = p
    )
  }
  if (predictive > 0 & !any_indep) {
    p <- base_dreamer_predictive_plot(
      x = x,
      doses = range_of_doses,
      time = times,
      probs = probs,
      reference_dose = reference_dose,
      predictive = predictive,
      p = p
    )
  } else if (predictive > 0 & any_indep) {
    p <- suppressMessages(
      base_dreamer_plot_error_bar(
        x = x,
        doses = doses,
        probs = probs,
        time = times,
        reference_dose = reference_dose,
        predictive = predictive,
        width = width * .5,
        p = p,
        color = "forestgreen",
        linetype = 2
      )
    )
  }
  p <- plot_data(p = p, data = data, x = x, times = times, doses = doses)
  if (!is.null(reference_dose)) {
    p <- p + labs(y = paste0("Response, relative to dose ", reference_dose))
  }
  return(p)
}

#' @title Plot Prior
#' @description Plot the prior over the dose range.  This is intended to
#'   help the user choose appropriate priors.
#' @param doses a vector of doses at which to evaluate and interpolate
#'   between.
#' @param probs A vector of length 2 indicating the lower and upper percentiles
#'   to plot.  Not applicable when `plot_draws = TRUE`.
#' @param n_samples the number of MCMC samples per MCMC chain used to generate
#'   the plot.
#' @param n_chains the number of MCMC chains.
#' @param ... model objects.  See \code{\link[dreamer]{model}} and
#'   examples below.
#' @param plot_draws if `TRUE`, the individual draws from the prior are plotted.
#'   If `FALSE`, only the prior mean and quantiles are drawn.
#' @param alpha the transparency setting for the prior draws in (0, 1].
#'   Only applies if `plot_draws = TRUE`.
#' @param times a vector of times at which to plot the prior.
#' @return The ggplot object.
#' @example man/examples/ex-dreamer_plot_prior.R
#' @export
dreamer_plot_prior <- function(
  n_samples = 1e4,
  probs = c(.025, .975),
  doses,
  n_chains = 1,
  ...,
  times = NULL,
  plot_draws = FALSE,
  alpha = .2
) {
  x <- dreamer_mcmc(
    data = NULL,
    n_adapt = 0,
    n_burn = 0,
    n_iter = n_samples,
    n_chains = n_chains,
    silent = TRUE,
    convergence_warn = FALSE,
    ...
  )
  any_longitudinal <- vapply(
      x,
      function(y) !is.null(attr(y, "longitudinal_model")),
      logical(1)
    ) %>%
    any()
  check_times(times, any_longitudinal)
  mods <- vapply(
    x,
    function(y) {
      any(inherits(y, c("dreamer_mcmc_continuous", "dreamer_mcmc_binary")))
    },
    logical(1)
  ) %>%
    which()
  for (i in mods) {
    attr(x[[i]], "doses") <- doses
  }
  any_independent <- vapply(
    x,
    function(y) inherits(y, "dreamer_mcmc_independent"),
    logical(1)
  ) %>%
    any()
  if (plot_draws) {
    names(doses) <- paste0("V", seq_len(length(doses)))
    samps <- posterior(
      x = x,
      doses = doses,
      times = times,
      probs = probs,
      return_samples = TRUE
    )$samps
    if (length(times) > 1) {
      samps <- samps %>%
        dplyr::mutate(
          dose = factor(.data$dose),
          iter_dose = paste0("iter_", .data$iter, "_dose", .data$dose)
        )
      p <- ggplot(
        samps,
        aes_string(
          x = "time",
          y = "mean_response",
          group = "iter_dose",
          color = "dose"
        )
      )
    } else {
      p <- ggplot(
        samps,
        aes_string(x = "dose", y = "mean_response", group = "iter")
      )
    }
    p <- p +
      geom_line(alpha = alpha) +
      labs(title = "Prior Draws", y = "Response")
  } else {
    p <- plot(x, doses = doses, times = times)
    if (length(times) == 1) {
      p <- p +
        labs(title = paste0(get_title(probs, "Prior"), ", time = ", times))
    }
  }
  if (length(times) > 1) {
    xlab <- "Time"
  } else {
    xlab <- "Dose"
  }
  p <- p + labs(
    x = xlab,
    y = "Response"
  ) +
    theme(plot.title = element_text(hjust = 0.5))
  return(p)
}

#' @title Compare Posterior Fits
#' @param ... `dreamer_mcmc` objects to be used for plotting.
#' @inheritParams dreamerplot
#' @param times the times at which to do the comparison.
#' @return a ggplot object.
#' @details If a Bayesian model averaging object is supplied first, all
#'   individual fits and the Bayesian model averaging fit will be plotted, with
#'   the model averaging fit in black (other model colors specified in the
#'   legend).  Otherwise,
#'   named arguments must be supplied for each model, and only the models
#'   provided will be plotted.
#' @example man/examples/ex-plot_comparison.R
#' @export
plot_comparison <- function(..., doses, times, probs, data, n_smooth, width) {
  UseMethod("plot_comparison", list(...)[[1]])
}

#' @rdname plot_comparison
#' @export
plot_comparison.default <- function(
  ...,
  doses = attr(list(...)[[1]], "doses"),
  times = NULL,
  probs = c(.025, .975),
  data = NULL,
  n_smooth = 50,
  width = bar_width(doses)
) {
  if (is.null(times) & !is.null(attr(list(...)[[1]], "times"))) {
    times <- max(attr(list(...)[[1]], "times"))
  }
  p <- plot_comparison_worker(
    x = list(...),
    doses = doses,
    times = times,
    probs = probs,
    data = data,
    n_smooth = n_smooth,
    width = width
  )
  return(p)
}

#' @rdname plot_comparison
#' @export
plot_comparison.dreamer_bma <- function(
  ...,
  doses = x$doses,
  times = NULL,
  probs = c(.025, .975),
  data = NULL,
  n_smooth = 50,
  width = bar_width(doses)
) {
  x <- list(...)[[1]]
  times <- get_time(x, times, max_length = Inf)
  model_index <- vapply(
    x,
    function(model) {
      any(inherits(model, c("dreamer_mcmc_continuous", "dreamer_mcmc_binary")))
    },
    logical(1)
  ) %>%
    which()
  p <- plot_comparison_worker(
    x = x[model_index],
    doses = doses,
    times = times,
    probs = probs,
    data = data,
    n_smooth = n_smooth,
    width = width
  )
  # add model averaging
  range_of_doses <- seq(
    from = min(doses),
    to = max(doses),
    length.out = n_smooth
  )
  dat <- posterior(
    x = x,
    doses = range_of_doses,
    times = times,
    probs = probs,
    return_samples = FALSE,
  )$stats
  dat$lb <- dat[[sprintf("%.2f%%", 100 * probs[1])]]
  dat$ub <- dat[[sprintf("%.2f%%", 100 * probs[2])]]
  dat_error_bar <- posterior(
    x = x,
    doses = doses,
    times = times,
    probs = probs,
    return_samples = FALSE,
  )$stats
  if (!is.null(times)) {
    dat_error_bar <- dat_error_bar %>%
      dplyr::filter(.data$time %in% !!times)
  }
  dat_error_bar$lb <- dat_error_bar[[sprintf("%.2f%%", 100 * probs[1])]]
  dat_error_bar$ub <- dat_error_bar[[sprintf("%.2f%%", 100 * probs[2])]]
  if (length(times) > 1) {
    xvar <- "time"
  } else {
    xvar <- "dose"
  }
  p <- p +
    geom_line(
      data = dat,
      mapping = aes_string(x = xvar, y = "mean"),
      inherit.aes = FALSE,
      lwd = 1.25
    ) +
    geom_line(
      data = dat,
      mapping = aes_string(x = xvar, y = "lb"),
      linetype = 2,
      inherit.aes = FALSE,
      lwd = 1.25
    ) +
    geom_line(
      data = dat,
      mapping = aes_string(x = xvar, y = "ub"),
      linetype = 2,
      inherit.aes = FALSE,
      lwd = 1.25
    ) +
    geom_errorbar(
      data = dat_error_bar,
      mapping = aes_string(x = xvar, ymin = "lb", ymax = "ub"),
      width = width,
      lwd = 1.25,
      inherit.aes = FALSE
    )
  return(p)
}

plot_comparison_worker <- function(
  x,
  doses,
  times,
  probs,
  data,
  n_smooth,
  width
) {
  check_names(x)
  dat_error_bar <- purrr::map(
    x,
    posterior,
    probs = probs,
    return_samples = FALSE,
    doses = doses,
    times = times
  ) %>%
    lapply(function(y) y$stats) %>%
    dplyr::bind_rows(.id = "Model")
  if (!is.null(times)) {
    dat_error_bar <- dat_error_bar %>%
      dplyr::filter(.data$time %in% !!times)
  }
  dat_error_bar$lb <- dat_error_bar[[sprintf("%.2f%%", 100 * probs[1])]]
  dat_error_bar$ub <- dat_error_bar[[sprintf("%.2f%%", 100 * probs[2])]]
  if (length(times) > 1) {
    assert_dose_len(doses)
    range_of_times <- seq(
      from = min(times),
      to = max(times),
      length.out = n_smooth
    )
    dat <- purrr::map(
      x,
      posterior,
      probs = probs,
      return_samples = FALSE,
      doses = doses,
      times = range_of_times
    ) %>%
      lapply(function(y) y$stats) %>%
      dplyr::bind_rows(.id = "Model")
    dat$lb <- dat[[sprintf("%.2f%%", 100 * probs[1])]]
    dat$ub <- dat[[sprintf("%.2f%%", 100 * probs[2])]]
    p <- ggplot(
      data = dat,
      mapping = aes_string(group = "Model", color = "Model")
    ) +
      geom_line(mapping = aes_string(x = "time", y = "mean")) +
      geom_line(mapping = aes_string(x = "time", y = "lb"), linetype = 2) +
      geom_line(mapping = aes_string(x = "time", y = "ub"), linetype = 2) +
      geom_errorbar(
        data = dat_error_bar,
        mapping = aes_string(x = "time", ymin = "lb", ymax = "ub"),
        width = width
      ) +
      labs(
        title = paste0(get_title(probs, "Posterior"), ", dose = ", doses),
        x = "Time",
        y = "Response"
      ) +
      theme(plot.title = element_text(hjust = 0.5))
  } else {
    range_of_doses <- seq(
      from = min(doses),
      to = max(doses),
      length.out = n_smooth
    )
    dat <- purrr::map(
      x,
      posterior,
      probs = probs,
      return_samples = FALSE,
      doses = range_of_doses,
      times = times
    ) %>%
      lapply(function(y) y$stats) %>%
      dplyr::bind_rows(.id = "Model")
    dat$lb <- dat[[sprintf("%.2f%%", 100 * probs[1])]]
    dat$ub <- dat[[sprintf("%.2f%%", 100 * probs[2])]]
    the_title <- get_title(probs, "Posterior") %>%
      add_time_to_title(times)
    p <- ggplot(
      data = dat,
      mapping = aes_string(group = "Model", color = "Model")
    ) +
      geom_line(mapping = aes_string(x = "dose", y = "mean")) +
      geom_line(mapping = aes_string(x = "dose", y = "lb"), linetype = 2) +
      geom_line(mapping = aes_string(x = "dose", y = "ub"), linetype = 2) +
      geom_errorbar(
        data = dat_error_bar,
        mapping = aes_string(x = "dose", ymin = "lb", ymax = "ub"),
        width = width
      ) +
      labs(
        title = the_title,
        x = "Dose",
        y = "Response"
      ) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  p <- plot_data(p = p, data = data, x = x[[1]], times = times, doses = doses)
  return(p)
}

plot_longitudinal <- function(
  x,
  doses,
  times,
  probs,
  n_smooth,
  predictive,
  reference_dose,
  data
) {
  width <- bar_width(times)
  post_bars <- posterior(
    x,
    doses = doses,
    times = times,
    reference_dose = reference_dose,
    probs = probs
  )$stats %>%
    dplyr::mutate(
      lb = .[[sprintf("%.2f%%", 100 * probs[1])]],
      ub = .[[sprintf("%.2f%%", 100 * probs[2])]],
      dose = factor(.data$dose)
    )
  times_smooth <- seq(from = min(times), to = max(times), length.out = n_smooth)
  post_smooth <- posterior(
    x,
    doses = doses,
    times = times_smooth,
    reference_dose = reference_dose,
    probs = probs
  )$stats %>%
    dplyr::mutate(
      lb = .[[sprintf("%.2f%%", 100 * probs[1])]],
      ub = .[[sprintf("%.2f%%", 100 * probs[2])]],
      dose = factor(.data$dose)
    )
  p <- ggplot(
    post_smooth,
    aes_string(x = "time", group = "dose", color = "dose")
  ) +
    geom_errorbar(
      data = post_bars,
      aes_string(ymin = "lb", ymax = "ub"),
      width = width
    ) +
    geom_line(aes_string(y = "mean")) +
    geom_line(aes_string(y = "lb"), lty = 2) +
    geom_line(aes_string(y = "ub"), lty = 2) +
    scale_x_continuous(breaks = times) +
    labs(x = "Time", y = "Response") +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(title = get_title(probs, "Posterior"))
  if (predictive > 0) {
    post_smooth_pred <- posterior(
      x,
      doses = doses,
      times = times_smooth,
      reference_dose = reference_dose,
      predictive = predictive,
      probs = probs
    )$stats %>%
      dplyr::mutate(
        lb = .[[sprintf("%.2f%%", 100 * probs[1])]],
        ub = .[[sprintf("%.2f%%", 100 * probs[2])]],
        dose = factor(.data$dose)
      )
    p <- p +
      geom_line(data = post_smooth_pred, aes_string(y = "lb"), lty = 3) +
      geom_line(data = post_smooth_pred, aes_string(y = "ub"), lty = 3)
  }
  p <- plot_data(p, data, x, times, doses)
  return(p)
}

check_plot <- function(doses, probs) {
  if (length(probs) != 2) {
    stop("Length(probs) must be 2", call. = FALSE)
  }
  if (is.null(doses)) {
    stop("Please provide doses", call. = FALSE)
  }
}

get_title <- function(probs, txt) {
  paste0(txt,
         " mean (solid) and ",
         sprintf("%.1f%%, %.1f%%", 100 * probs[1], 100 * probs[2]),
         " quantiles (dashed)"
  )
}

get_title_error_bars <- function(probs) {
  paste0(
    sprintf("%.1f%%, %.1f%%", 100 * probs[1], 100 * probs[2]),
    " Posterior Quantiles"
  )
}

bar_width <- function(doses) {
  if (length(doses) > 1) {
    width <- diff(range(doses)) / 15
  } else {
    width <- 1 / 10
  }
  return(width)
}

any_independent <- function(x) {
  UseMethod("any_independent", x)
}

any_independent.dreamer_bma <- function(x) {
  vapply(
    x,
    function(y) inherits(y, "dreamer_mcmc_independent"),
    logical(1)
  ) %>%
    any()
}

any_independent.default <- function(x) {
  inherits(x, "dreamer_mcmc_independent")
}

aggregate_data <- function(data, type) {
  if (type == "binary") {
    agdat <- dplyr::group_by(data, across(any_of(c("dose", "time")))) %>%
      dplyr::summarize(
        response = sum(.data$response) / sum(.data$n)
      )
  } else {
    agdat <- dplyr::group_by(data, across(any_of(c("dose", "time")))) %>%
      dplyr::summarize(
        response = sum(.data$response * .data$n) / sum(.data$n)
      )
  }
  return(agdat)
}

plot_data <- function(p, data, x, times, doses) {
  if (is.null(data)) return(p)
  if (!rlang::has_name(data, "n")) {
    group_vars <- "dose"
    if (!is.null(attr(x, "longitudinal_model"))) {
      group_vars <- c(group_vars, "time")
    }
    agdat <- data %>%
      dplyr::group_by(across(any_of(!!group_vars))) %>%
      dplyr::summarise(response = mean(.data$response)) %>%
      dplyr::ungroup()
  } else {
    agdat <- aggregate_data(
      data = data,
      type = attr(x, "response_type")
    )
  }
  if (rlang::has_name(agdat, "time")) {
    agdat <- agdat %>%
      dplyr::filter(.data$time %in% !!times)
  }
  agdat <- agdat %>%
    dplyr::filter(.data$dose %in% !!doses)
  if (length(times) > 1) {
    agdat <- agdat %>%
      dplyr::mutate(dose = factor(.data$dose, levels = levels(p$data$dose)))
    p <- p + geom_point(
      aes_string(x = "time", y = "response", color = "dose"),
      data = agdat,
      size = 2,
      inherit.aes = FALSE
    )
  } else {
    p <- p + geom_point(
      aes_string(x = "dose", y = "response"),
      data = agdat,
      size = 2,
      inherit.aes = FALSE
    )
  }
  return(p)
}

check_no_dots <- function(function_name, ...) {
  vars <- list(...)
  if (length(vars) > 0) {
    stop(
      paste0(
        "Extra (unused) arguments in ",
        function_name,
        ": ",
        paste0(names(vars), collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

check_times <- function(times, any_longitudinal) {
  if (any_longitudinal & is.null(times)) {
    stop(
      "Must specify 'times' argument with longitudinal models.",
      call. = FALSE
    )
  }
  if (!any_longitudinal & !is.null(times)) {
    stop(
      "Can only specify 'times' argument with longitudinal models.",
      call. = FALSE
    )
  }
}

check_names <- function(x) {
  if (is.null(names(x)) | any(names(x) == "")) {
    stop("All models must be named", call. = FALSE)
  }
}

assert_dose_len <- function(doses) {
  if (length(doses) > 1) {
    stop("If length(times) > 1, length(doses) must be 1", call. = FALSE)
  }
}
