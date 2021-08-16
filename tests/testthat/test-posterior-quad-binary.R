test_that("MCMC: quad binary logit", {
  n_chains <- 2
  link <- "logit"
  data <- dreamer_data_linear_binary(
    n_cohorts = c(10, 20, 30),
    dose = c(1, 3, 5),
    b1 = 1,
    b2 = 2,
    link = link
  )
  mcmc <- dreamer_mcmc(
    data,
    mod = model_quad_binary(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      mu_b3 = 0,
      sigma_b3 = 1,
      link = link
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
    b1 = 1:10 / 100,
    b2 = 2:11 / 100,
    b3 = 3:12 / 100,
    true_responses = rlang::expr(ilogit(b1 + b2 * dose + b3 * dose ^ 2))
  )
  # with dose adjustment
  test_posterior(
    mcmc,
    doses = c(1, 3, 2),
    reference_dose = .5,
    prob = c(.25, .75),
    b1 = 1:10 / 100,
    b2 = 2:11 / 100,
    b3 = 3:12 / 100,
    true_responses = rlang::expr(
      ilogit(b1 + b2 * dose + b3 * dose ^ 2) -
        ilogit(b1 + b2 * reference_dose + b3 * reference_dose ^ 2)
    )
  )
})

test_that("MCMC: quad binary probit", {
  n_chains <- 2
  link <- "probit"
  data <- dreamer_data_linear_binary(
    n_cohorts = c(10, 20, 30),
    dose = c(1, 3, 5),
    b1 = 1,
    b2 = 2,
    link = link
  )
  mcmc <- dreamer_mcmc(
    data,
    mod = model_quad_binary(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      mu_b3 = 0,
      sigma_b3 = 1,
      link = link
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
    b1 = 1:10 / 100,
    b2 = 2:11 / 100,
    b3 = 3:12 / 100,
    true_responses = rlang::expr(iprobit(b1 + b2 * dose + b3 * dose ^ 2))
  )
  # with dose adjustment
  test_posterior(
    mcmc,
    doses = c(1, 3, 2),
    reference_dose = .5,
    prob = c(.25, .75),
    b1 = 1:10 / 100,
    b2 = 2:11 / 100,
    b3 = 3:12 / 100,
    true_responses = rlang::expr(
      iprobit(b1 + b2 * dose + b3 * dose ^ 2) -
        iprobit(b1 + b2 * reference_dose + b3 * reference_dose ^ 2)
    )
  )
})

test_that("MCMC: quad binary logit long linear", {
  n_chains <- 2
  t_max <- 4
  times <- c(0, 2, 4)
  link <- "logit"
  data <- dreamer_data_linear_binary(
    n_cohorts = c(10, 20, 30),
    dose = c(1, 3, 5),
    b1 = 1,
    b2 = 2,
    link = link,
    longitudinal = "linear",
    a = .5,
    times = times,
    t_max = t_max
  )
  mcmc <- dreamer_mcmc(
    data,
    mod = model_quad_binary(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      mu_b3 = 0,
      sigma_b3 = 1,
      link = link,
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
    b1 = 1:10 / 100,
    b2 = 2:11 / 100,
    b3 = 3:12 / 100,
    a = 10:1 / 100,
    true_responses = rlang::expr(
      ilogit(a + time / !!t_max * (b1 + b2 * dose + b3 * dose ^ 2))
    )
  )
  test_posterior(
    mcmc,
    doses = c(1, 3, 2),
    times = c(1, 5, 2),
    reference_dose = .5,
    prob = c(.25, .75),
    b1 = 1:10 / 100,
    b2 = 2:11 / 100,
    b3 = 3:12 / 100,
    a = 10:1 / 100,
    true_responses = rlang::expr(
      ilogit(a + (time / !!t_max) * (b1 + b2 * dose + b3 * dose ^ 2)) -
        ilogit(
          a +
            time / !!t_max *
            (b1 + b2 * reference_dose + b3 * reference_dose ^ 2)
        )
    )
  )
})

test_that("MCMC: quad binary probit long linear", {
  n_chains <- 2
  t_max <- 4
  times <- c(0, 2, 4)
  link <- "probit"
  data <- dreamer_data_linear_binary(
    n_cohorts = c(10, 20, 30),
    dose = c(1, 3, 5),
    b1 = 1,
    b2 = 2,
    link = link,
    longitudinal = "linear",
    a = .5,
    times = times,
    t_max = t_max
  )
  mcmc <- dreamer_mcmc(
    data,
    mod = model_quad_binary(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      mu_b3 = 0,
      sigma_b3 = 1,
      link = link,
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
    b1 = 1:10 / 100,
    b2 = 2:11 / 100,
    b3 = 3:12 / 100,
    a = 10:1 / 100,
    true_responses = rlang::expr(
      iprobit(a + time / !!t_max * (b1 + b2 * dose + b3 * dose ^ 2))
    )
  )
  test_posterior(
    mcmc,
    doses = c(1, 3, 2),
    times = c(1, 5, 2),
    reference_dose = .5,
    prob = c(.25, .75),
    b1 = 1:10 / 100,
    b2 = 2:11 / 100,
    b3 = 3:12 / 100,
    a = 10:1 / 100,
    true_responses = rlang::expr(
      iprobit(a + (time / !!t_max) * (b1 + b2 * dose + b3 * dose ^ 2)) -
        iprobit(
          a +
            time / !!t_max *
            (b1 + b2 * reference_dose + b3 * reference_dose ^ 2)
        )
    )
  )
})


test_that("MCMC: quad binary logit long ITP", {
  n_chains <- 2
  t_max <- 4
  times <- c(0, 2, 4)
  link <- "logit"
  data <- dreamer_data_linear_binary(
    n_cohorts = c(10, 20, 30),
    dose = c(1, 3, 5),
    b1 = 1,
    b2 = 2,
    link = link,
    longitudinal = "linear",
    a = .5,
    times = times,
    t_max = t_max
  )
  mcmc <- dreamer_mcmc(
    data,
    mod = model_quad_binary(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      mu_b3 = 0,
      sigma_b3 = 1,
      link = link,
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
    b1 = 1:10 / 100,
    b2 = 2:11 / 100,
    b3 = 3:12 / 100,
    a = 10:1 / 100,
    c1 = seq(.1, 3, length = 10) / 100,
    true_responses = rlang::expr(
      ilogit(
        a + (1 - exp(- c1 * time)) / (1 - exp(- c1 * !!t_max)) *
          (b1 + b2 * dose + b3 * dose ^ 2)
      )
    )
  )
  test_posterior(
    mcmc,
    doses = c(1, 3, 2),
    times = c(1, 5, 2),
    reference_dose = .5,
    prob = c(.25, .75),
    b1 = 1:10 / 100,
    b2 = 2:11 / 100,
    b3 = 3:12 / 100,
    a = 10:1 / 100,
    c1 = seq(.1, 3, length = 10) / 100,
    true_responses = rlang::expr(
      ilogit(
        a + (1 - exp(- c1 * time)) / (1 - exp(- c1 * !!t_max)) *
          (b1 + b2 * dose + b3 * dose ^ 2)
      ) -
        ilogit(
          (a + (1 - exp(- c1 * time)) / (1 - exp(- c1 * !!t_max)) *
             (b1 + b2 * reference_dose + b3 * reference_dose ^ 2))
        )
    )
  )
})

test_that("MCMC: quad binary probit long ITP", {
  n_chains <- 2
  t_max <- 4
  times <- c(0, 2, 4)
  link <- "probit"
  data <- dreamer_data_linear_binary(
    n_cohorts = c(10, 20, 30),
    dose = c(1, 3, 5),
    b1 = 1,
    b2 = 2,
    link = link,
    longitudinal = "linear",
    a = .5,
    times = times,
    t_max = t_max
  )
  mcmc <- dreamer_mcmc(
    data,
    mod = model_quad_binary(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      mu_b3 = 0,
      sigma_b3 = 1,
      link = link,
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
    b1 = 1:10 / 100,
    b2 = 2:11 / 100,
    b3 = 3:12 / 100,
    a = 10:1 / 100,
    c1 = seq(.1, 3, length = 10) / 100,
    true_responses = rlang::expr(
      iprobit(
        a + (1 - exp(- c1 * time)) / (1 - exp(- c1 * !!t_max)) *
          (b1 + b2 * dose + b3 * dose ^ 2)
      )
    )
  )
  test_posterior(
    mcmc,
    doses = c(1, 3, 2),
    times = c(1, 5, 2),
    reference_dose = .5,
    prob = c(.25, .75),
    b1 = 1:10 / 100,
    b2 = 2:11 / 100,
    b3 = 3:12 / 100,
    a = 10:1 / 100,
    c1 = seq(.1, 3, length = 10) / 100,
    true_responses = rlang::expr(
      iprobit(
        a + (1 - exp(- c1 * time)) / (1 - exp(- c1 * !!t_max)) *
          (b1 + b2 * dose + b3 * dose ^ 2)
      ) -
        iprobit(
          (a + (1 - exp(- c1 * time)) / (1 - exp(- c1 * !!t_max)) *
             (b1 + b2 * reference_dose + b3 * reference_dose ^ 2))
        )
    )
  )
})

test_that("MCMC: quad binary logit long IDP", {
  n_chains <- 2
  t_max <- 4
  times <- c(0, 2, 4)
  link <- "logit"
  data <- dreamer_data_linear_binary(
    n_cohorts = c(10, 20, 30),
    dose = c(1, 3, 5),
    b1 = 1,
    b2 = 2,
    link = link,
    longitudinal = "linear",
    a = .5,
    times = times,
    t_max = t_max
  )
  mcmc <- dreamer_mcmc(
    data,
    mod = model_quad_binary(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      mu_b3 = 0,
      sigma_b3 = 1,
      link = link,
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
    b1 = 1:10 / 100,
    b2 = 2:11 / 100,
    b3 = 3:12 / 100,
    a = 10:1 / 100,
    c1 = seq(.1, 3, length = 10) / 100,
    c2 = seq(-.1, -.02, length = 10) / 100,
    d1 = seq(3, 4, length = 10) / 100,
    d2 = seq(4, 5, length = 10) / 100,
    gam = seq(.2, .33, length = 10) / 100,
    true_responses = rlang::expr(
      ilogit(
        a + (b1 + b2 * dose + b3 * dose ^ 2) * (
          (1 - exp(- c1 * time)) / (1 - exp(- c1 * d1)) * (time < d1) +
            (1 - gam * (1 - exp(- c2 * (time - d1))) /
               (1 - exp(- c2 * (d2 - d1)))) * (d1 <= time & time <= d2) +
            (1 - gam) * (time > d2)
        )
      )
    )
  )
  test_posterior(
    mcmc,
    doses = c(1, 3, 2),
    times = c(1, 5, 2),
    reference_dose = .5,
    prob = c(.25, .75),
    a = 10:1 / 100,
    b1 = 1:10 / 100,
    b2 = 2:11 / 100,
    b3 = 3:12 / 100,
    c1 = seq(.1, 3, length = 10) / 100,
    c2 = seq(-.1, -.02, length = 10) / 100,
    d1 = seq(3, 4, length = 10) / 100,
    d2 = seq(4, 5, length = 10) / 100,
    gam = seq(.2, .33, length = 10) / 100,
    true_responses = rlang::expr(
      ilogit(
        a + (b1 + b2 * dose + b3 * dose ^ 2) * (
          (1 - exp(- c1 * time)) / (1 - exp(- c1 * d1)) * (time < d1) +
            (1 - gam * (1 - exp(- c2 * (time - d1))) /
               (1 - exp(- c2 * (d2 - d1)))) * (d1 <= time & time <= d2) +
            (1 - gam) * (time > d2)
        )
      ) -
        ilogit(
          (
            a + (b1 + b2 * reference_dose + b3 * reference_dose ^ 2) * (
              (1 - exp(- c1 * time)) / (1 - exp(- c1 * d1)) * (time < d1) +
                (1 - gam * (1 - exp(- c2 * (time - d1))) /
                   (1 - exp(- c2 * (d2 - d1)))) * (d1 <= time & time <= d2) +
                (1 - gam) * (time > d2)
            )
          )
        )
    )
  )
})

test_that("MCMC: quad binary probit long IDP", {
  n_chains <- 2
  t_max <- 4
  times <- c(0, 2, 4)
  link <- "probit"
  data <- dreamer_data_linear_binary(
    n_cohorts = c(10, 20, 30),
    dose = c(1, 3, 5),
    b1 = 1,
    b2 = 2,
    link = link,
    longitudinal = "linear",
    a = .5,
    times = times,
    t_max = t_max
  )
  mcmc <- dreamer_mcmc(
    data,
    mod = model_quad_binary(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      mu_b3 = 0,
      sigma_b3 = 1,
      link = link,
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
    b1 = 1:10 / 100,
    b2 = 2:11 / 100,
    b3 = 3:12 / 100,
    a = 10:1 / 100,
    c1 = seq(.1, 3, length = 10) / 100,
    c2 = seq(-.1, -.02, length = 10) / 100,
    d1 = seq(3, 4, length = 10) / 100,
    d2 = seq(4, 5, length = 10) / 100,
    gam = seq(.2, .33, length = 10) / 100,
    true_responses = rlang::expr(
      iprobit(
        a + (b1 + b2 * dose + b3 * dose ^ 2) * (
          (1 - exp(- c1 * time)) / (1 - exp(- c1 * d1)) * (time < d1) +
            (1 - gam * (1 - exp(- c2 * (time - d1))) /
               (1 - exp(- c2 * (d2 - d1)))) * (d1 <= time & time <= d2) +
            (1 - gam) * (time > d2)
        )
      )
    )
  )
  test_posterior(
    mcmc,
    doses = c(1, 3, 2),
    times = c(1, 5, 2),
    reference_dose = .5,
    prob = c(.25, .75),
    a = 10:1 / 100,
    b1 = 1:10 / 100,
    b2 = 2:11 / 100,
    b3 = 3:12 / 100,
    c1 = seq(.1, 3, length = 10) / 100,
    c2 = seq(-.1, -.02, length = 10) / 100,
    d1 = seq(3, 4, length = 10) / 100,
    d2 = seq(4, 5, length = 10) / 100,
    gam = seq(.2, .33, length = 10) / 100,
    true_responses = rlang::expr(
      iprobit(
        a + (b1 + b2 * dose + b3 * dose ^ 2) * (
          (1 - exp(- c1 * time)) / (1 - exp(- c1 * d1)) * (time < d1) +
            (1 - gam * (1 - exp(- c2 * (time - d1))) /
               (1 - exp(- c2 * (d2 - d1)))) * (d1 <= time & time <= d2) +
            (1 - gam) * (time > d2)
        )
      ) -
        iprobit(
          (
            a + (b1 + b2 * reference_dose + b3 * reference_dose ^ 2) * (
              (1 - exp(- c1 * time)) / (1 - exp(- c1 * d1)) * (time < d1) +
                (1 - gam * (1 - exp(- c2 * (time - d1))) /
                   (1 - exp(- c2 * (d2 - d1)))) * (d1 <= time & time <= d2) +
                (1 - gam) * (time > d2)
            )
          )
        )
    )
  )
})
