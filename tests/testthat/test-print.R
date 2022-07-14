test_that("print methods", {
  set.seed(12)
  data <- dreamer_data_linear(
    n_cohorts = c(10, 20, 30),
    dose = c(0, 2.5, 5),
    b1 = 1,
    b2 = 2,
    sigma = 3
  )
  mod <- model_linear(
    mu_b1 = 0,
    sigma_b1 = 1,
    mu_b2 = 0,
    sigma_b2 = 1,
    shape = 1,
    rate = .001,
    w_prior = 1
  )
  output <- dreamer_mcmc(
    data = data,
    n_adapt = 1e3,
    n_burn = 0,
    n_iter = 10,
    n_chains = 1,
    silent = TRUE,
    convergence_warn = FALSE,
    mod_linear = mod
  )
  expect_snapshot(print(mod))
  expect_snapshot(print(output))
  expect_snapshot(print(output$mod_linear))
  
  mod_binary <-  model_linear_binary(
    mu_b1 = 0,
    sigma_b1 = 1,
    mu_b2 = 0,
    sigma_b2 = 1,
    link = "probit",
    w_prior = 1
  )
  expect_snapshot(print(mod_binary))
})
