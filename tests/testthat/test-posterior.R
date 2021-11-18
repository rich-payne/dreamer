test_that("posterior.dreamer_bma() model_index and model names", {
  set.seed(885)
  t_max <- 5
  data <- dreamer_data_linear(
    n_cohorts = c(10, 15, 20, 25, 30),
    dose = c(1:5),
    b1 = .5,
    b2 = 3,
    sigma = .5,
    times = 1:5,
    longitudinal = "linear",
    a = 5,
    t_max = t_max
  )
  out <- dreamer_mcmc(
    data = data,
    n_iter = 10,
    n_chains = 2,
    convergence_warn = FALSE,
    silent = TRUE,
    mod_lin = model_linear(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      shape = 1,
      rate = .01,
      w_prior = .5
    ),
    mod_quad = model_quad(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      mu_b3 = 0,
      sigma_b3 = 1,
      shape = 1,
      rate = .01,
      w_prior = .5
    )
  )
  model_index <- attr(out, "model_index")
  samps_lin <- posterior(out$mod_lin, return_samples = TRUE)$samps
  samps_quad <- posterior(out$mod_quad, return_samples = TRUE)$samps
  samps <- dplyr::bind_rows(
    samps_lin[model_index == 1, ],
    samps_quad[model_index == 2, ]
  ) %>%
    dplyr::arrange(iter, dose)
  obs <- posterior(out, return_samples = TRUE)$samps %>%
    dplyr::arrange(iter, dose)
  expect_equal(obs, samps)
})

test_that("posterior.dreamer_bma uses model_index and iter correctly", {
  set.seed(885)
  t_max <- 5
  data <- dreamer_data_linear(
    n_cohorts = c(10, 15, 20, 25, 30),
    dose = c(1:5),
    b1 = .5,
    b2 = 3,
    sigma = .5,
    times = 1:5,
    longitudinal = "linear",
    a = 5,
    t_max = t_max
  )
  out <- dreamer_mcmc(
    data = data,
    n_iter = 10,
    n_chains = 2,
    convergence_warn = FALSE,
    silent = TRUE,
    mod_lin = model_linear(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      shape = 1,
      rate = .01,
      w_prior = .5
    ),
    mod_quad = model_quad(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      mu_b3 = 0,
      sigma_b3 = 1,
      shape = 1,
      rate = .01,
      w_prior = .5
    )
  )
  iter <- c(1, 3:10)
  model_index <- attr(out, "model_index")
  samps_lin <- posterior(out$mod_lin, return_samples = TRUE)$samps
  samps_quad <- posterior(out$mod_quad, return_samples = TRUE)$samps
  samps <- dplyr::bind_rows(
    samps_lin[model_index == 1, ],
    samps_quad[model_index == 2, ]
  ) %>%
    dplyr::arrange(iter, dose) %>%
    dplyr::filter(iter %in% !!iter)
  obs <- posterior(out, return_samples = TRUE, iter = iter)$samps %>%
    dplyr::arrange(iter, dose)
  expect_equal(obs, samps)
})

test_that("continuous predictive runs", {
  data <- dreamer_data_linear(n_cohorts = c(10, 20, 30), c(1, 3, 5), 1, 2, 2)
  mcmc <- dreamer_mcmc(
    data,
    mod = model_linear(0, 1, 0, 1, 1, 1),
    n_iter = 10,
    silent = TRUE,
    convergence_warn = FALSE,
    n_chains = 1
  )
  expect_true(
    tibble::is_tibble(
      posterior(mcmc, reference_dose = 0, predictive = 10)$stats
    )
  )
})

test_that("binary predictive runs", {
  data <- dreamer_data_linear_binary(
    n_cohorts = c(10, 20, 30),
    dose = c(1, 3, 5),
    b1 = 1,
    b2 = 2,
    link = "logit"
  )
  mcmc <- dreamer_mcmc(
    data,
    mod = model_linear_binary(0, 1, 0, 1, link = "logit"),
    n_iter = 10,
    silent = TRUE,
    convergence_warn = FALSE,
    n_chains = 1
  )
  stats <- posterior(mcmc, reference_dose = 0, predictive = 10)$stats
  expect_true(tibble::is_tibble(stats))
  expect_equal(
    colnames(stats),
    c("dose", "reference_dose", "mean", "2.50%", "97.50%")
  )
  stats_pred <- posterior(mcmc, predictive = 10)$stats
  expect_true(tibble::is_tibble(stats_pred))
  expect_equal(
    colnames(stats_pred),
    c("dose", "mean", "2.50%", "97.50%")
  )
})

test_that("predictive works correctly (binary)", {
  set.seed(0)
  data <- dreamer_data_linear_binary(
    n_cohorts = c(10, 20, 30),
    dose = c(1, 3, 5),
    b1 = -1,
    b2 = .2,
    link = "logit"
  )
  mcmc <- dreamer_mcmc(
    data,
    mod = model_linear_binary(0, 1, 0, 1, link = "logit"),
    n_iter = 10,
    silent = TRUE,
    convergence_warn = FALSE,
    n_chains = 1
  )
  draw_pred_samps <- function(mcmc, predictive, rep) {
    samps <- posterior(
      mcmc,
      predictive = predictive,
      return_samples = TRUE,
      return_stats = FALSE
    )$samps$mean_response %>% mean()
  }
  predictive <- 5
  reps <- purrr::map_dbl(
    1:1000,
    draw_pred_samps,
    mcmc = mcmc,
    predictive = predictive
  )
  obs_var <- var(reps)
  exp_var <- posterior(
    mcmc,
    return_samples = TRUE,
    return_stats = FALSE
  )$samps %>%
    dplyr::summarize(
      true_var = sum((mean_response * (1 - mean_response) / !!predictive)) /
        n() ^ 2
    ) %>%
    dplyr::pull(true_var)

  # bootstrap to get range
  bstrap_qtiles <- quantile(
    replicate(1e3, var(sample(reps, length(reps), replace = TRUE))),
    c(.025, .975)
  )
  expect_true(bstrap_qtiles[1] < exp_var && exp_var < bstrap_qtiles[2])

  n_reps <- 1e3
  reps_xbar <- replicate(n_reps, posterior(mcmc, predictive = 5)$stats$mean)
  avg_reps_xbar <- apply(reps_xbar, 1, mean)
  se_reps_xbar <- apply(reps_xbar, 1, sd) / sqrt(n_reps)
  xbar_diff <- posterior(mcmc)$stats$mean - avg_reps_xbar
  expect_true(all(abs(xbar_diff) < 1.96 * se_reps_xbar))
})

test_that("predictive works correctly (binary) with dose adjustment", {
  set.seed(100)
  data <- dreamer_data_linear_binary(
    n_cohorts = c(10, 20, 30),
    dose = c(1, 3, 5),
    b1 = -1,
    b2 = .2,
    link = "logit"
  )
  mcmc <- dreamer_mcmc(
    data,
    mod = model_linear_binary(0, 1, 0, 1, link = "logit"),
    n_iter = 10,
    silent = TRUE,
    convergence_warn = FALSE,
    n_chains = 1
  )
  reference_dose <- 1.5
  draw_pred_samps <- function(mcmc, predictive, reference_dose, rep) {
    posterior(
      mcmc,
      predictive = predictive,
      reference_dose = reference_dose,
      return_samples = TRUE,
      return_stats = FALSE
    )$samps$mean_response %>% mean()
  }
  predictive <- 5
  reps <- purrr::map_dbl(
    1:1000,
    draw_pred_samps,
    mcmc = mcmc,
    predictive = predictive,
    reference_dose = reference_dose
  )
  obs_var <- var(reps)
  exp_var <- posterior(
    mcmc,
    return_samples = TRUE,
    return_stats = FALSE
  )$samps
  exp_var_ref <- posterior(
    mcmc,
    doses = reference_dose,
    return_samples = TRUE,
    return_stats = FALSE
  )$samps
  exp_var <- dplyr::left_join(
    exp_var,
    exp_var_ref,
    by = c("iter", "model"),
    suffix = c("", "_ref")
  )
  exp_var <- exp_var %>%
    dplyr::mutate(
      var1 = mean_response * (1 - mean_response) / (!!predictive * n() ^ 2),
      var2 = mean_response_ref * (1 - mean_response_ref) /
        (!!predictive * n() ^ 2),
      true_var = var1 + var2
    ) %>%
    dplyr::pull(true_var) %>%
    sum()
  # bootstrap to get range
  bstrap_qtiles <- quantile(
    replicate(1e3, var(sample(reps, length(reps), replace = TRUE))),
    c(.025, .975)
  )
  expect_true(bstrap_qtiles[1] < exp_var && exp_var < bstrap_qtiles[2])

  n_reps <- 1e3
  reps_xbar <- replicate(
    n_reps,
    posterior(
      mcmc,
      predictive = predictive,
      reference_dose = reference_dose
    )$stats$mean
  )
  avg_reps_xbar <- apply(reps_xbar, 1, mean)
  se_reps_xbar <- apply(reps_xbar, 1, sd) / sqrt(n_reps)
  xbar_diff <- posterior(mcmc, reference_dose = reference_dose)$stats$mean -
    avg_reps_xbar
  expect_true(all(abs(xbar_diff) < 1.96 * se_reps_xbar))
})

test_that("predictive works correctly (continuous)", {
  set.seed(10)
  data <- dreamer_data_linear(
    n_cohorts = c(10, 20, 30),
    dose = c(1, 3, 5),
    b1 = -1,
    b2 = .2,
    sigma = .25
  )
  mcmc <- dreamer_mcmc(
    data,
    mod = model_linear(0, 1, 0, 1, shape = 1, rate = .001),
    n_iter = 10,
    silent = TRUE,
    convergence_warn = FALSE,
    n_chains = 1
  )
  draw_pred_samps <- function(mcmc, predictive, rep) {
    samps <- posterior(
      mcmc,
      dose = 1,
      predictive = predictive,
      return_samples = TRUE,
      return_stats = FALSE
    )$samps$mean_response %>% mean()
  }
  predictive <- 5
  reps <- purrr::map_dbl(
    1:1000,
    draw_pred_samps,
    mcmc = mcmc,
    predictive = predictive
  )
  obs_var <- var(reps)
  exp_var <- tibble::as_tibble(as.matrix(mcmc$mod)) %>%
    dplyr::summarize(
      true_var = (sum(sigma ^ 2) / !!predictive) / n() ^ 2
    ) %>%
    dplyr::pull(true_var)

  # bootstrap to get range
  bstrap_qtiles <- quantile(
    replicate(1e3, var(sample(reps, length(reps), replace = TRUE))),
    c(.025, .975)
  )
  expect_true(bstrap_qtiles[1] < exp_var && exp_var < bstrap_qtiles[2])

  # ensure mean is the same
  n_reps <- 1e3
  reps_xbar <- replicate(
    n_reps,
    posterior(mcmc, predictive = predictive)$stats$mean
  )
  avg_reps_xbar <- apply(reps_xbar, 1, mean)
  se_reps_xbar <- apply(reps_xbar, 1, sd) / sqrt(n_reps)
  xbar_diff <- posterior(mcmc)$stats$mean - avg_reps_xbar
  expect_true(all(abs(xbar_diff) < 1.96 * se_reps_xbar))
})

test_that("predictive works correctly (continuous) adjusted", {
  set.seed(25)
  data <- dreamer_data_linear(
    n_cohorts = c(10, 20, 30),
    dose = c(1, 3, 5),
    b1 = -1,
    b2 = .2,
    sigma = .25
  )
  mcmc <- dreamer_mcmc(
    data,
    mod = model_linear(0, 1, 0, 1, shape = 1, rate = .001),
    n_iter = 10,
    silent = TRUE,
    convergence_warn = FALSE,
    n_chains = 1
  )
  draw_pred_samps <- function(mcmc, predictive, rep, reference_dose = .5) {
    samps <- posterior(
      mcmc,
      dose = 1,
      reference_dose = reference_dose,
      predictive = predictive,
      return_samples = TRUE,
      return_stats = FALSE
    )$samps$mean_response %>% mean()
  }
  predictive <- 5
  reps <- purrr::map_dbl(
    1:1000,
    draw_pred_samps,
    mcmc = mcmc,
    predictive = predictive
  )
  obs_var <- var(reps)
  exp_var <- tibble::as_tibble(as.matrix(mcmc$mod)) %>%
    dplyr::summarize(
      true_var = (sum(sigma ^ 2) / (!!predictive / 2)) / n() ^ 2
    ) %>%
    dplyr::pull(true_var)

  # bootstrap to get range
  bstrap_qtiles <- quantile(
    replicate(1e3, var(sample(reps, length(reps), replace = TRUE))),
    c(.025, .975)
  )
  expect_true(bstrap_qtiles[1] < exp_var && exp_var < bstrap_qtiles[2])

  # ensure mean is the same
  n_reps <- 1e3
  reps_xbar <- replicate(
    n_reps,
    posterior(mcmc, predictive = predictive)$stats$mean
  )
  avg_reps_xbar <- apply(reps_xbar, 1, mean)
  se_reps_xbar <- apply(reps_xbar, 1, sd) / sqrt(n_reps)
  xbar_diff <- posterior(mcmc)$stats$mean - avg_reps_xbar
  expect_true(all(abs(xbar_diff) < 1.96 * se_reps_xbar))
})
