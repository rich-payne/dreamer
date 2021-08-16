test_that("get_extreme_monotonic()", {
  data <- dreamer_data_linear(n_cohorts = c(10, 20, 30), c(1, 3, 5), 1, 2, 2)
  mcmc <- dreamer_mcmc(
    data,
    mod = model_linear(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
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
  b2 <- c(-1, 1)
  mcmc <- mcmc %>%
    replace_mcmc("mod", "b1", b1) %>%
    replace_mcmc("mod", "b2", b2)
  obs <- get_extreme_monotonic(
    mcmc$mod,
    time = NULL,
    greater = TRUE,
    lower = lower,
    upper = upper,
    index = NULL
  )
  exp <- tibble::tibble(doses = c(lower, upper)) %>%
    dplyr::mutate(
      extreme_responses = b1 + b2 * doses,
      greater = TRUE
    )
  expect_equal(obs, exp)

  obs <- get_extreme(
    mcmc$mod,
    time = NULL,
    greater = FALSE,
    lower = lower,
    upper = upper,
    index = NULL
  )
  exp <- tibble::tibble(doses = c(upper, lower)) %>%
    dplyr::mutate(
      extreme_responses = b1 + b2 * doses,
      greater = FALSE
    )
  expect_equal(obs, exp)

  obs <- get_extreme(
    mcmc$mod,
    time = NULL,
    greater = FALSE,
    lower = lower,
    upper = upper,
    index = 2
  )
  exp <- tibble::tibble(doses = c(upper, lower)) %>%
    dplyr::mutate(
      extreme_responses = b1 + b2 * doses,
      greater = FALSE
    ) %>%
    dplyr::slice(2)
  expect_equal(obs, exp)
})

test_that("get_extreme_monotonic() longitudinal", {
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
    mod = model_linear(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      shape = 1,
      rate = .01,
      longitudinal = model_longitudinal_linear(0, 1, t_max)
    )
  )
  lower <- min(attr(mcmc, "doses"))
  upper <- max(attr(mcmc, "doses"))
  a <- c(.1, .2)
  b1 <- c(1, 5)
  b2 <- c(-1, 1)

  mcmc <- mcmc %>%
    replace_mcmc("mod", "a", a) %>%
    replace_mcmc("mod", "b1", b1) %>%
    replace_mcmc("mod", "b2", b2)

  time <- 3
  obs <- get_extreme_monotonic(
    mcmc$mod,
    time = time,
    greater = TRUE,
    lower = lower,
    upper = upper,
    index = NULL
  )
  exp <- tibble::tibble(doses = c(lower, upper)) %>%
    dplyr::mutate(
      extreme_responses = a + (b1 + b2 * doses) * time / t_max,
      greater = TRUE
    )
  expect_equal(obs, exp)

  obs <- get_extreme(
    mcmc$mod,
    time = time,
    greater = FALSE,
    lower = lower,
    upper = upper,
    index = NULL
  )
  exp <- tibble::tibble(doses = c(upper, lower)) %>%
    dplyr::mutate(
      extreme_responses = a + (b1 + b2 * doses) * time / t_max,
      greater = FALSE
    )
  expect_equal(obs, exp)

  obs <- get_extreme(
    mcmc$mod,
    time = time,
    greater = FALSE,
    lower = lower,
    upper = upper,
    index = 2
  )
  exp <- tibble::tibble(doses = c(upper, lower)) %>%
    dplyr::mutate(
      extreme_responses = a + (b1 + b2 * doses) * time / t_max,
      greater = FALSE
    ) %>%
    dplyr::slice(2)
  expect_equal(obs, exp)
})
