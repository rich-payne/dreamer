test_that("global seed implies same MCMC", {
  set.seed(888)
  data <- dreamer_data_quad(
    n_cohorts = c(10, 10, 10, 10),
    dose = c(.25, .5, .75, 1.5),
    b1 = 0,
    b2 = 2,
    b3 = -1,
    sigma = .5
  )

  set.seed(1)
  output1 <- dreamer_mcmc(
    data = data,
    n_adapt = 1e3,
    n_burn = 1e3,
    n_iter = 1e3,
    n_chains = 2,
    silent = TRUE,
    mod_linear = model_linear(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      shape = 1,
      rate = .001,
      w_prior = 1
    )
  )

  set.seed(1)
  output2 <- dreamer_mcmc(
    data = data,
    n_adapt = 1e3,
    n_burn = 1e3,
    n_iter = 1e3,
    n_chains = 2,
    silent = TRUE,
    mod_linear = model_linear(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      shape = 1,
      rate = .001,
      w_prior = 1
    )
  )

  expect_identical(output1, output2)
})

test_that("different seeds imply different MCMC", {
  set.seed(888)
  data <- dreamer_data_quad(
    n_cohorts = c(10, 10, 10, 10),
    dose = c(.25, .5, .75, 1.5),
    b1 = 0,
    b2 = 2,
    b3 = -1,
    sigma = .5
  )

  set.seed(1)
  output1 <- dreamer_mcmc(
    data = data,
    n_adapt = 1e3,
    n_burn = 1e3,
    n_iter = 1e3,
    n_chains = 2,
    silent = TRUE,
    mod_linear = model_linear(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      shape = 1,
      rate = .001,
      w_prior = 1
    )
  )

  set.seed(2)
  output2 <- dreamer_mcmc(
    data = data,
    n_adapt = 1e3,
    n_burn = 1e3,
    n_iter = 1e3,
    n_chains = 2,
    silent = TRUE,
    mod_linear = model_linear(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      shape = 1,
      rate = .001,
      w_prior = 1
    )
  )

  expect_false(identical(output1, output2))
})
