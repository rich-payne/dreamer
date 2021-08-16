set.seed(88552211)
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
  n_iter = 100,
  n_chains = 2,
  silent = TRUE,
  convergence_warn = FALSE
)

test_that("summary() runs", {
  out <- summary(mcmc)
  expect_true(tibble::is_tibble(out$model_weights))
  expect_true(tibble::is_tibble(out$summary))
})

test_that("diagnostics() runs", {
  out <- diagnostics(mcmc)
  expect_true(tibble::is_tibble(out))
})

test_that("summary() with one parameter", {
  data <- dreamer_data_linear_binary(n_cohorts = 10, 1, 1, 2, link = "logit")
  mcmc <- dreamer_mcmc(
    data,
    mod = model_independent_binary(
      mu_b1 = 0,
      sigma_b1 = 1,
      link = "logit"
    ),
    n_iter = 100,
    n_chains = 2,
    silent = TRUE,
    convergence_warn = FALSE
  )
  out <- summary(mcmc)
  expect_true(tibble::is_tibble(out$model_weights))
  expect_true(tibble::is_tibble(out$summary))
})
