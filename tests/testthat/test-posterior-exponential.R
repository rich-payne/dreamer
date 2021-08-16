test_that("MCMC: exponential", {
  n_chains <- 2
  data <- dreamer_data_linear(n_cohorts = c(10, 20, 30), c(1, 3, 5), 1, 2, 2)
  mcmc <- dreamer_mcmc(
    data,
    mod = model_exp(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      mu_b3 = 0,
      sigma_b3 = 1,
      shape = 1,
      rate = .01
    ),
    n_iter = 5,
    silent = TRUE,
    convergence_warn = FALSE,
    n_chains = n_chains
  )
  assert_mcmc_format(mcmc, n_chains)
  # dreamer post
  test_posterior(
    mcmc,
    doses = c(1, 3, 2),
    prob = c(.25, .75),
    b1 = 1:10,
    b2 = 2:11,
    b3 = 3:12,
    true_responses = rlang::expr(b1 + b2 * (1 - exp(- b3 * dose)))
  )
  # with dose adjustment
  test_posterior(
    mcmc,
    doses = c(1, 3, 2),
    reference_dose = .5,
    prob = c(.25, .75),
    b1 = 1:10,
    b2 = 2:11,
    b3 = 3:12,
    true_responses = rlang::expr(
      b1 + b2 * (1 - exp(- b3 * dose)) -
        (b1 + b2 * (1 - exp(- b3 * reference_dose)))
    )
  )
})

test_that("MCMC: exponential long linear", {
  n_chains <- 2
  t_max <- 4
  times <- c(0, 2, 4)
  data <- dreamer_data_linear(
    n_cohorts = c(10, 20, 30),
    dose = c(1, 3, 5),
    b1 = 1,
    b2 = 2,
    sigma = 2,
    longitudinal = "linear",
    a = .5,
    times = times,
    t_max = t_max
  )
  mcmc <- dreamer_mcmc(
    data,
    mod = model_exp(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      mu_b3 = 0,
      sigma_b3 = 1,
      shape = 1,
      rate = .01,
      longitudinal = model_longitudinal_linear(0, 1, t_max)
    ),
    n_iter = 5,
    silent = TRUE,
    convergence_warn = FALSE,
    n_chains = n_chains
  )

  assert_mcmc_format(mcmc, n_chains, times)
  test_posterior(
    mcmc,
    doses = c(1, 3, 2),
    times = c(1, 5, 2),
    prob = c(.25, .75),
    b1 = 1:10,
    b2 = 2:11,
    b3 = 3:12,
    a = 10:1,
    true_responses = rlang::expr(
      a + time / !!t_max * (b1 + b2 * (1 - exp(- b3 * dose)))
    )
  )
  test_posterior(
    mcmc,
    doses = c(1, 3, 2),
    times = c(1, 5, 2),
    reference_dose = .5,
    prob = c(.25, .75),
    b1 = 1:10,
    b2 = 2:11,
    b3 = 3:12,
    a = 10:1,
    true_responses = rlang::expr(
      a + (time / !!t_max) * (b1 + b2 * (1 - exp(- b3 * dose))) -
        (a + (time / !!t_max) * (b1 + b2 * (1 - exp(- b3 * reference_dose))))
    )
  )
})

test_that("MCMC: exponential long ITP", {
  n_chains <- 2
  t_max <- 4
  times <- c(0, 2, 4)
  data <- dreamer_data_linear(
    n_cohorts = c(10, 20, 30),
    dose = c(1, 3, 5),
    b1 = 1,
    b2 = 2,
    sigma = 2,
    longitudinal = "linear",
    a = .5,
    times = times,
    t_max = t_max
  )
  mcmc <- dreamer_mcmc(
    data,
    mod = model_exp(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      mu_b3 = 0,
      sigma_b3 = 1,
      shape = 1,
      rate = .01,
      longitudinal = model_longitudinal_itp(0, 1, t_max = t_max)
    ),
    n_iter = 5,
    silent = TRUE,
    convergence_warn = FALSE,
    n_chains = n_chains
  )

  assert_mcmc_format(mcmc, n_chains, times)
  test_posterior(
    mcmc,
    doses = c(1, 3, 2),
    times = c(1, 5, 2),
    prob = c(.25, .75),
    b1 = 1:10,
    b2 = 2:11,
    b3 = 3:12,
    a = 10:1,
    c1 = seq(.1, 3, length = 10),
    true_responses = rlang::expr(
      a + (1 - exp(- c1 * time)) / (1 - exp(- c1 * !!t_max)) *
        (b1 + b2 * (1 - exp(- b3 * dose)))
    )
  )
  test_posterior(
    mcmc,
    doses = c(1, 3, 2),
    times = c(1, 5, 2),
    reference_dose = .5,
    prob = c(.25, .75),
    b1 = 1:10,
    b2 = 2:11,
    b3 = 3:12,
    a = 10:1,
    c1 = seq(.1, 3, length = 10),
    true_responses = rlang::expr(
      a + (1 - exp(- c1 * time)) / (1 - exp(- c1 * !!t_max)) *
        (b1 + b2 * (1 - exp(- b3 * dose))) -
        (a + (1 - exp(- c1 * time)) / (1 - exp(- c1 * !!t_max)) *
           (b1 + b2 * (1 - exp(- b3 * reference_dose))))
    )
  )
})

test_that("MCMC: exponential long IDP", {
  n_chains <- 2
  t_max <- 4
  times <- c(0, 2, 4)
  data <- dreamer_data_linear(
    n_cohorts = c(10, 20, 30),
    dose = c(1, 3, 5),
    b1 = 1,
    b2 = 2,
    sigma = 2,
    longitudinal = "linear",
    a = .5,
    times = times,
    t_max = t_max
  )
  mcmc <- dreamer_mcmc(
    data,
    mod = model_exp(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      mu_b3 = 0,
      sigma_b3 = 1,
      shape = 1,
      rate = .01,
      longitudinal = model_longitudinal_idp(0, 1, t_max = t_max)
    ),
    n_iter = 5,
    silent = TRUE,
    convergence_warn = FALSE,
    n_chains = n_chains
  )

  assert_mcmc_format(mcmc, n_chains, times)
  test_posterior(
    mcmc,
    doses = c(1, 3, 2),
    times = c(1, 5, 2),
    prob = c(.25, .75),
    b1 = 1:10,
    b2 = 2:11,
    b3 = 3:12,
    a = 10:1,
    c1 = seq(.1, 3, length = 10),
    c2 = seq(-.1, -.02, length = 10),
    d1 = seq(3, 4, length = 10),
    d2 = seq(4, 5, length = 10),
    gam = seq(.2, .33, length = 10),
    true_responses = rlang::expr(
      a + (b1 + b2 * (1 - exp(- b3 * dose))) * (
        (1 - exp(- c1 * time)) / (1 - exp(- c1 * d1)) * (time < d1) +
          (1 - gam * (1 - exp(- c2 * (time - d1))) /
             (1 - exp(- c2 * (d2 - d1)))) * (d1 <= time & time <= d2) +
          (1 - gam) * (time > d2)
      )
    )
  )
  test_posterior(
    mcmc,
    doses = c(1, 3, 2),
    times = c(1, 5, 2),
    reference_dose = .5,
    prob = c(.25, .75),
    a = 10:1,
    b1 = 1:10,
    b2 = 2:11,
    b3 = 3:12,
    c1 = seq(.1, 3, length = 10),
    c2 = seq(-.1, -.02, length = 10),
    d1 = seq(3, 4, length = 10),
    d2 = seq(4, 5, length = 10),
    gam = seq(.2, .33, length = 10),
    true_responses = rlang::expr(
      a + (b1 + b2 * (1 - exp(- b3 * dose))) * (
        (1 - exp(- c1 * time)) / (1 - exp(- c1 * d1)) * (time < d1) +
          (1 - gam * (1 - exp(- c2 * (time - d1))) /
             (1 - exp(- c2 * (d2 - d1)))) * (d1 <= time & time <= d2) +
          (1 - gam) * (time > d2)
      ) -
        (
          a + (b1 + b2 * (1 - exp(- b3 * reference_dose))) * (
            (1 - exp(- c1 * time)) / (1 - exp(- c1 * d1)) * (time < d1) +
              (1 - gam * (1 - exp(- c2 * (time - d1))) /
                 (1 - exp(- c2 * (d2 - d1)))) * (d1 <= time & time <= d2) +
              (1 - gam) * (time > d2)
          )
        )
    )
  )
})
