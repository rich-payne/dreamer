test_that("get_dose() beta", {
  data <- dreamer_data_linear(n_cohorts = c(10, 20, 30), c(1, 3, 5), 1, 2, 2)
  mcmc <- dreamer_mcmc(
    data,
    mod = model_beta(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      mu_b3 = 0,
      sigma_b3 = 1,
      mu_b4 = 0,
      sigma_b4 = 1,
      shape = 1,
      rate = .01
    ),
    n_iter = 2,
    n_chains = 1,
    silent = TRUE,
    convergence_warn = FALSE
  )
  lower <- min(attr(mcmc, "doses"))
  upper <- max(attr(mcmc, "doses"))

  b1 <- 1:2
  b2 <- c(- 1, 1)
  b3 <- c(2.2, 2)
  b4 <- c(.99, 1.01)

  mcmc <- mcmc %>%
    replace_mcmc("mod", "b1", b1) %>%
    replace_mcmc("mod", "b2", b2) %>%
    replace_mcmc("mod", "b3", b3) %>%
    replace_mcmc("mod", "b4", b4)

  scale <- attr(mcmc$mod, "scale")
  dose <- 2
  get_dose(
    mcmc$mod,
    time = NULL,
    response = (b1 + b2 * ((b3 + b4) ^ (b3 + b4)) / (b3 ^ b3 * b4 ^ b4) *
                  (dose / scale) ^ b3 * (1 - dose / scale) ^ b4),
    lower = lower,
    upper = upper
  ) %>%
    expect_equal(rep(dose, 2))
})

test_that("get_dose() beta longitudinal", {
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
    data = data,
    n_iter = 2,
    n_chains = 1,
    convergence_warn = FALSE,
    silent = TRUE,
    mod = model_beta(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      mu_b3 = 0,
      sigma_b3 = 1,
      mu_b4 = 0,
      sigma_b4 = 1,
      shape = 1,
      rate = .01,
      longitudinal = model_longitudinal_linear(0, 1, t_max)
    )
  )
  lower <- min(attr(mcmc, "doses"))
  upper <- max(attr(mcmc, "doses"))

  a <- c(.1, .2)
  b1 <- 1:2
  b2 <- c(- 1, 1)
  b3 <- c(1.98, 2)
  b4 <- c(.99, 1.01)

  mcmc <- mcmc %>%
    replace_mcmc("mod", "a", a) %>%
    replace_mcmc("mod", "b1", b1) %>%
    replace_mcmc("mod", "b2", b2) %>%
    replace_mcmc("mod", "b3", b3) %>%
    replace_mcmc("mod", "b4", b4)

  time <- 3
  dose <- 2
  scale <- attr(mcmc$mod, "scale")
  get_dose(
    mcmc$mod,
    time = time,
    response =
      a + time / t_max *
        ((b1 + b2 * ((b3 + b4) ^ (b3 + b4)) / (b3 ^ b3 * b4 ^ b4) *
            (dose / scale) ^ b3 * (1 - dose / scale) ^ b4)),
    lower = lower,
    upper = upper
  ) %>%
    expect_equal(rep(dose, 2))
})
