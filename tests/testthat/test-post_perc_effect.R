test_that("post_perc_effect.dreamer()", {
  data <- dreamer_data_linear(n_cohorts = c(10, 20, 30), c(1, 3, 5), 1, 2, 2)
  mcmc <- dreamer_mcmc(
    data,
    mod = model_linear(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      shape = 1,
      rate = .01
    ),
    n_iter = 10,
    silent = TRUE,
    convergence_warn = FALSE,
    n_chains = 1
  )

  small_bound <- .05
  lower <- 0
  upper <- 5
  b1 <- 1:10
  b2 <- c(-c(5:1), 1:5)
  mcmc <- mcmc %>%
    replace_mcmc("mod", "b1", b1) %>%
    replace_mcmc("mod", "b2", b2)

  dose <- 2.5
  max_dose_gt <- c(rep(lower, 5), rep(upper, 5))
  ed100_gt <- b1 + b2 * max_dose_gt
  post_percs_gt <- ((b1 + b2 * dose) - small_bound) / (ed100_gt - small_bound)
  post_percs_gt[!(post_percs_gt >= 0 & post_percs_gt <= 1)] <- NA

  exp <- tibble::tibble(
    dose,
    pr_perc_exists = mean(!is.na(post_percs_gt)),
    mean = mean(post_percs_gt, na.rm = TRUE),
    `2.5%` = quantile(post_percs_gt, prob = .025, na.rm = TRUE),
    `97.5%` = quantile(post_percs_gt, prob = .975, na.rm = TRUE)
  )
  obs <- post_perc_effect(
    mcmc,
    dose,
    small_bound = small_bound,
    lower = lower,
    upper = upper
  )$stats
  expect_equal(obs, exp)

  max_dose_lt <- c(rep(upper, 5), rep(lower, 5))
  ed100_lt <- b1 + b2 * max_dose_lt
  post_percs_lt <- ((b1 + b2 * dose) - small_bound) / (ed100_lt - small_bound)
  post_percs_lt[!(post_percs_lt >= 0 & post_percs_lt <= 1)] <- NA

  exp <- tibble::tibble(
    dose,
    pr_perc_exists = mean(!is.na(post_percs_lt)),
    mean = mean(post_percs_lt, na.rm = TRUE),
    `2.5%` = quantile(post_percs_lt, prob = .025, na.rm = TRUE),
    `97.5%` = quantile(post_percs_lt, prob = .975, na.rm = TRUE)
  )
  obs <- post_perc_effect(
    mcmc,
    dose,
    greater = FALSE,
    small_bound = small_bound,
    lower = lower,
    upper = upper
  )$stats
  expect_equal(obs, exp)
})

test_that("post_perc_effect.dreamer() longitudinal", {
  times <- c(0, 10)
  t_max <- max(times)
  data <- dreamer_data_linear(
    n_cohorts = c(10, 25, 30),
    dose = c(0, 2, 4),
    b1 = .5,
    b2 = 3,
    sigma = .5,
    longitudinal = "linear",
    a = .5,
    times = times,
    t_max = t_max
  )
  mcmc <- dreamer_mcmc(
    data,
    mod = model_linear(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      shape = 1,
      rate = .01,
      longitudinal = model_longitudinal_linear(0, 1, t_max)
    ),
    n_iter = 10,
    silent = TRUE,
    convergence_warn = FALSE,
    n_chains = 1
  )

  small_bound <- .05
  lower <- 0
  upper <- 5
  a <- c(0:9) / 10
  b1 <- 1:10
  b2 <- c(-c(5:1), 1:5)
  mcmc <- mcmc %>%
    replace_mcmc("mod", "a", a) %>%
    replace_mcmc("mod", "b1", b1) %>%
    replace_mcmc("mod", "b2", b2)

  dose <- 2.5
  time <- 4
  max_dose_gt <- c(rep(lower, 5), rep(upper, 5))
  ed100_gt <- a + time / t_max * (b1 + b2 * max_dose_gt)
  post_percs_gt <- (a + time / t_max * (b1 + b2 * dose) - small_bound) /
    (ed100_gt - small_bound)
  post_percs_gt[!(post_percs_gt >= 0 & post_percs_gt <= 1)] <- NA

  exp <- tibble::tibble(
    dose,
    pr_perc_exists = mean(!is.na(post_percs_gt)),
    mean = mean(post_percs_gt, na.rm = TRUE),
    `2.5%` = quantile(post_percs_gt, prob = .025, na.rm = TRUE),
    `97.5%` = quantile(post_percs_gt, prob = .975, na.rm = TRUE)
  )
  obs <- post_perc_effect(
    mcmc,
    dose,
    time = time,
    small_bound = small_bound,
    lower = lower,
    upper = upper
  )$stats
  expect_equal(obs, exp)

  max_dose_lt <- c(rep(upper, 5), rep(lower, 5))
  ed100_lt <- a + time / t_max * (b1 + b2 * max_dose_lt)
  post_percs_lt <- (a + time / t_max * (b1 + b2 * dose) - small_bound) /
    (ed100_lt - small_bound)
  post_percs_lt[!(post_percs_lt >= 0 & post_percs_lt <= 1)] <- NA

  exp <- tibble::tibble(
    dose,
    pr_perc_exists = mean(!is.na(post_percs_lt)),
    mean = mean(post_percs_lt, na.rm = TRUE),
    `2.5%` = quantile(post_percs_lt, prob = .025, na.rm = TRUE),
    `97.5%` = quantile(post_percs_lt, prob = .975, na.rm = TRUE)
  )
  obs <- post_perc_effect(
    mcmc,
    dose,
    time = time,
    greater = FALSE,
    small_bound = small_bound,
    lower = lower,
    upper = upper
  )$stats
  expect_equal(obs, exp)
})

test_that("post_perc_effect.dreamer_bma()", {
  data <- dreamer_data_linear(n_cohorts = c(10, 20, 30), c(1, 3, 5), 1, 2, 2)
  set.seed(88332)
  mcmc <- dreamer_mcmc(
    data,
    mod = model_linear(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      shape = 1,
      rate = .01,
      w_prior = .5
    ),
    quad = model_quad(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      mu_b3 = 0,
      sigma_b3 = 1,
      shape = 1,
      rate = .01,
      w_prior = .5
    ),
    n_iter = 5,
    silent = TRUE,
    convergence_warn = FALSE,
    n_chains = 1
  )

  model_index <- attr(mcmc, "model_index")
  dose <- 3.5
  exp <- dplyr::bind_rows(
    post_perc_effect(mcmc$mod, dose = dose, return_samples = TRUE)$samps %>%
      dplyr::slice(which(model_index == 1)),
    post_perc_effect(mcmc$quad, dose = dose, return_samples = TRUE)$samps %>%
      dplyr::slice(which(model_index == 2))
  ) %>%
    dplyr::group_by(dose) %>%
    dplyr::summarize(
      pr_perc_exists = mean(!is.na(post_eds)),
      mean = mean(post_eds),
      `2.5%` = quantile(post_eds, prob = .025),
      `97.5%` = quantile(post_eds, prob = .975),
      .groups = "drop"
    )
  expect_equal(post_perc_effect(mcmc, dose = dose)$stats, exp)
})
