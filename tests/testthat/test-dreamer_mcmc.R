test_that("model_index is sampled correctly", {
  data <- dreamer_data_linear(
    n_cohorts = c(10, 20, 30),
    dose = c(0, 2.5, 5),
    b1 = 1,
    b2 = 2,
    sigma = 3
  )

  output <- dreamer_mcmc(
    data = data,
    n_adapt = 1e3,
    n_burn = 1e3,
    n_iter = 1e3,
    n_chains = 1,
    silent = TRUE,
    convergence_warn = FALSE,
    mod_linear = model_linear(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      shape = 1,
      rate = .001,
      w_prior = 0
    ),
    mod_logquad = model_logquad(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      mu_b3 = 0,
      sigma_b3 = 1,
      shape = 1,
      rate = .001,
      w_prior = 1
    )
  )
  expect_true(all(attr(output, "model_index") == 2))
})

test_that("model names attributes", {
  data <- dreamer_data_linear(
    n_cohorts = c(10, 20, 30),
    dose = c(0, 2.5, 5),
    b1 = 1,
    b2 = 2,
    sigma = 3
  )
  output <- dreamer_mcmc(
    data = data,
    n_adapt = 1e3,
    n_burn = 1e3,
    n_iter = 1e3,
    n_chains = 1,
    silent = TRUE,
    convergence_warn = FALSE,
    mod_logquad = model_logquad(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      mu_b3 = 0,
      sigma_b3 = 1,
      shape = 1,
      rate = .001,
      w_prior = .5
    ),
    mod_linear = model_linear(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      shape = 1,
      rate = .001,
      w_prior = .5
    )
  )
  expect_equal(attr(output, "model_names"), c("mod_logquad", "mod_linear"))
  expect_equal(attr(output$mod_linear, "model_name"), "mod_linear")
  expect_equal(attr(output$mod_logquad, "model_name"), "mod_logquad")
})

test_that("jags modules are restored after MCMC", {
  data <- dreamer_data_linear(
    n_cohorts = c(10, 20, 30),
    dose = c(0, 2.5, 5),
    b1 = 1,
    b2 = 2,
    sigma = 3
  )
  original_modules <- rjags::list.modules()
  output <- dreamer_mcmc(
    data = data,
    n_adapt = 1e3,
    n_burn = 1e3,
    n_iter = 1e3,
    n_chains = 1,
    silent = TRUE,
    convergence_warn = FALSE,
    mod_linear = model_linear(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      shape = 1,
      rate = .001,
      w_prior = 0
    ),
    mod_logquad = model_logquad(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      mu_b3 = 0,
      sigma_b3 = 1,
      shape = 1,
      rate = .001,
      w_prior = 1
    )
  )
  expect_equal(rjags::list.modules(), original_modules)
})
