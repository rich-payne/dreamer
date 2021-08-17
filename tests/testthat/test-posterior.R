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
    n_iter = 5,
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
    n_iter = 5,
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
    n_iter = 5,
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
    n_iter = 5,
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
