test_that("dreamer_mean_longitudinal_mcmc", {
  set.seed(3226262)
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
      longitudinal = model_longitudinal_idp(0, 1, 0, 1, 0, 1, t_max = t_max)
    ),
    n_iter = 10,
    silent = TRUE,
    convergence_warn = FALSE,
    n_chains = 1
  )

  mcmc_mat <- as.matrix(mcmc$mod)
  samps <- tibble::as_tibble(mcmc_mat)

  time <- 7.3
  dreamer_mean_longitudinal_mcmc(
    time = time,
    longitudinal_model = NULL,
    mcmc_mat,
    t_max
  ) %>%
    expect_equal(1)

  dreamer_mean_longitudinal_mcmc(
    time = time,
    longitudinal_model = "dreamer_longitudinal_linear",
    mcmc_mat,
    t_max
  ) %>%
    expect_equal(
      time / t_max
    )

  dreamer_mean_longitudinal_mcmc(
    time = time,
    longitudinal_model = "dreamer_longitudinal_itp",
    mcmc_mat,
    t_max
  ) %>%
    expect_equal(
      (1 - exp(- samps$c1 * time)) / (1 - exp(- samps$c1 * t_max))
    )

  dreamer_mean_longitudinal_mcmc(
    time = time,
    longitudinal_model = "dreamer_longitudinal_idp",
    mcmc_mat,
    t_max
  ) %>%
    expect_equal(
      (1 - exp(- samps$c1 * time)) / (1 - exp(- samps$c1 * samps$d1)) *
        (time < samps$d1) +
          (1 - samps$gam * (1 - exp(- samps$c2 * (time - samps$d1))) /
             (1 - exp(- samps$c2 * (samps$d2 - samps$d1)))) *
             (samps$d1 <= time & time <= samps$d2) +
        (1 - samps$gam) * (time > samps$d2)
    )
})
