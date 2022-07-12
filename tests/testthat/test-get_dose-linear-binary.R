test_that("get_dose() linear binary logit", {
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
    mod = model_linear_binary(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
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
  b2 <- c(- 1, 1)

  mcmc <- mcmc %>%
    replace_mcmc("mod", "b1", b1) %>%
    replace_mcmc("mod", "b2", b2)

  dose <- 2
  get_dose(
    mcmc$mod,
    time = NULL,
    response = ilogit(b1 + b2 * dose),
    lower = lower,
    upper = upper
  ) %>%
    expect_equal(rep(dose, 2))
})

test_that("get_dose() linear binary probit", {
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
    mod = model_linear_binary(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
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
  b2 <- c(- 1, 1)

  mcmc <- mcmc %>%
    replace_mcmc("mod", "b1", b1) %>%
    replace_mcmc("mod", "b2", b2)

  dose <- 2
  get_dose(
    mcmc$mod,
    time = NULL,
    response = iprobit(b1 + b2 * dose),
    lower = lower,
    upper = upper
  ) %>%
    expect_equal(rep(dose, 2))
})

test_that("get_dose() linear binary logit longitudinal", {
  times <- c(0, 10)
  t_max <- max(times)
  link <- "logit"
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
    mod = model_linear_binary(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      link = link,
      longitudinal = model_longitudinal_linear(0, 1, t_max)
    )
  )
  lower <- min(attr(mcmc, "doses"))
  upper <- max(attr(mcmc, "doses"))

  a <- c(.1, .2)
  b1 <- 1:2
  b2 <- c(- 1, 1)

  mcmc <- mcmc %>%
    replace_mcmc("mod", "a", a) %>%
    replace_mcmc("mod", "b1", b1) %>%
    replace_mcmc("mod", "b2", b2)

  time <- 3
  dose <- 2
  get_dose(
    mcmc$mod,
    time = time,
    response = ilogit(.2 + time / t_max * (2 + 1 * dose)),
    lower = lower,
    upper = upper
  ) %>%
    expect_equal(c(NA, dose))
})

test_that("get_dose() linear binary probit longitudinal", {
  times <- c(0, 10)
  t_max <- max(times)
  link <- "probit"
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
    mod = model_linear_binary(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      link = link,
      longitudinal = model_longitudinal_linear(0, 1, t_max)
    )
  )
  lower <- min(attr(mcmc, "doses"))
  upper <- max(attr(mcmc, "doses"))

  a <- c(.1, .2)
  b1 <- 1:2
  b2 <- c(- 1, 1)

  mcmc <- mcmc %>%
    replace_mcmc("mod", "a", a) %>%
    replace_mcmc("mod", "b1", b1) %>%
    replace_mcmc("mod", "b2", b2)

  time <- 3
  dose <- 2
  get_dose(
    mcmc$mod,
    time = time,
    response = iprobit(.2 + time / t_max * (2 + 1 * dose)),
    lower = lower,
    upper = upper
  ) %>%
    expect_equal(c(NA, dose))
})
