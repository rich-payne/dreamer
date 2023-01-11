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
  mod_binary <-  model_linear_binary(
    mu_b1 = 0,
    sigma_b1 = 1,
    mu_b2 = 0,
    sigma_b2 = 1,
    link = "probit",
    w_prior = 1
  )
  expect_output(print(mod))
  expect_output(print(output))
  expect_output(print(output$mod_linear))
  expect_output(print(mod_binary))
  skip_on_ci()
  skip_on_cran()
  expect_snapshot(print(mod))
  expect_snapshot(print(output))
  expect_snapshot(print(output$mod_linear))
  expect_snapshot(print(mod_binary))
})

test_that("print methods (longitudinal)", {
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
  mod <- model_linear(
    mu_b1 = 0,
    sigma_b1 = 1,
    mu_b2 = 0,
    sigma_b2 = 1,
    shape = 1,
    rate = .01,
    w_prior = 1,
    longitudinal = model_longitudinal_linear(0, 1, t_max)
  )
  out <- dreamer_mcmc(
    data = data,
    n_iter = 10,
    n_chains = 2,
    convergence_warn = FALSE,
    silent = TRUE,
    mod_lin = mod
  )
  expect_output(print(mod))
  expect_output(print(out))
  expect_output(print(out$mod_lin))
  skip_on_ci()
  skip_on_cran()
  expect_snapshot(print(mod))
  expect_snapshot(print(out))
  expect_snapshot(print(out$mod_lin))
})
