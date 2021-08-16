test_that("get_dose() emax binary logit", {
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
    mod = model_emax_binary(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      mu_b3 = 0,
      sigma_b3 = 1,
      mu_b4 = 0,
      sigma_b4 = 1,
      link = link
    ),
    n_iter = 2,
    n_chains = 1,
    silent = TRUE,
    convergence_warn = FALSE
  )
  lower <- min(attr(mcmc, "doses"))
  upper <- max(attr(mcmc, "doses"))

  b1 <- 1:2
  b2 <- c(-1, 1)
  b3 <- c(-2, 2)
  b4 <- c(.99, 1.01)

  mcmc <- mcmc %>%
    replace_mcmc("mod", "b1", b1) %>%
    replace_mcmc("mod", "b2", b2) %>%
    replace_mcmc("mod", "b3", b3) %>%
    replace_mcmc("mod", "b4", b4)

  dose <- 2
  get_dose(
    mcmc$mod,
    time = NULL,
    response = ilogit(
      b1 + (b2 - b1) * dose ^ b4 / (exp(b3 * b4) + dose ^ b4)
    ),
    lower = lower,
    upper = upper
  ) %>%
    expect_equal(rep(dose, 2))
})

test_that("get_dose() emax binary probit", {
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
    mod = model_emax_binary(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      mu_b3 = 0,
      sigma_b3 = 1,
      mu_b4 = 0,
      sigma_b4 = 1,
      link = link
    ),
    n_iter = 2,
    n_chains = 1,
    silent = TRUE,
    convergence_warn = FALSE
  )
  lower <- min(attr(mcmc, "doses"))
  upper <- max(attr(mcmc, "doses"))

  b1 <- 1:2
  b2 <- c(-1, 1)
  b3 <- c(-2, 2)
  b4 <- c(.99, 1.01)

  mcmc <- mcmc %>%
    replace_mcmc("mod", "b1", b1) %>%
    replace_mcmc("mod", "b2", b2) %>%
    replace_mcmc("mod", "b3", b3) %>%
    replace_mcmc("mod", "b4", b4)

  dose <- 2
  get_dose(
    mcmc$mod,
    time = NULL,
    response = iprobit(
      b1 + (b2 - b1) * dose ^ b4 / (exp(b3 * b4) + dose ^ b4)
    ),
    lower = lower,
    upper = upper
  ) %>%
    expect_equal(rep(dose, 2))
})

test_that("get_dose() emax binary logit longitudinal", {
  link <- "logit"
  times <- c(0, 10)
  t_max <- max(times)
  data <- dreamer_data_linear_binary(
    n_cohorts = c(10, 25, 30),
    dose = c(0, 2, 4),
    b1 = .5,
    b2 = 3,
    link = link,
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
    mod = model_emax_binary(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      mu_b3 = 0,
      sigma_b3 = 1,
      mu_b4 = 0,
      sigma_b4 = 1,
      link = link,
      longitudinal = model_longitudinal_linear(0, 1, t_max)
    )
  )
  lower <- min(attr(mcmc, "doses"))
  upper <- max(attr(mcmc, "doses"))

  a <- c(.1, .2)
  b1 <- 1:2
  b2 <- c(-1, 1)
  b3 <- c(-2, 2)
  b4 <- c(.99, 1.01)

  mcmc <- mcmc %>%
    replace_mcmc("mod", "a", a) %>%
    replace_mcmc("mod", "b1", b1) %>%
    replace_mcmc("mod", "b2", b2) %>%
    replace_mcmc("mod", "b3", b3) %>%
    replace_mcmc("mod", "b4", b4)

  time <- 3
  dose <- 2
  get_dose(
    mcmc$mod,
    time = time,
    response =
      ilogit(
        a + time / t_max *
        (b1 + (b2 - b1) * dose ^ b4 / (exp(b3 * b4) + dose ^ b4))
      ),
    lower = lower,
    upper = upper
  ) %>%
    expect_equal(rep(dose, 2))
})

test_that("get_dose() emax binary probit longitudinal", {
  link <- "probit"
  times <- c(0, 10)
  t_max <- max(times)
  data <- dreamer_data_linear_binary(
    n_cohorts = c(10, 25, 30),
    dose = c(0, 2, 4),
    b1 = .5,
    b2 = 3,
    link = link,
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
    mod = model_emax_binary(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      mu_b3 = 0,
      sigma_b3 = 1,
      mu_b4 = 0,
      sigma_b4 = 1,
      link = link,
      longitudinal = model_longitudinal_linear(0, 1, t_max)
    )
  )
  lower <- min(attr(mcmc, "doses"))
  upper <- max(attr(mcmc, "doses"))

  a <- c(.1, .2)
  b1 <- 1:2
  b2 <- c(-1, 1)
  b3 <- c(-2, 2)
  b4 <- c(.99, 1.01)

  mcmc <- mcmc %>%
    replace_mcmc("mod", "a", a) %>%
    replace_mcmc("mod", "b1", b1) %>%
    replace_mcmc("mod", "b2", b2) %>%
    replace_mcmc("mod", "b3", b3) %>%
    replace_mcmc("mod", "b4", b4)

  time <- 3
  dose <- 2
  get_dose(
    mcmc$mod,
    time = time,
    response =
      iprobit(
        a + time / t_max *
          (b1 + (b2 - b1) * dose ^ b4 / (exp(b3 * b4) + dose ^ b4))
      ),
    lower = lower,
    upper = upper
  ) %>%
    expect_equal(rep(dose, 2))
})
