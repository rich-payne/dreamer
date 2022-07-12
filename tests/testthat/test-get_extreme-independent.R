test_that("get_extreme.independent()", {
  doses <- c(1, 3, 5)
  data <- dreamer_data_linear(n_cohorts = c(10, 20, 30), doses, 1, 2, 2)
  mcmc <- dreamer_mcmc(
    data,
    mod = model_independent(
      mu_b1 = 0,
      sigma_b1 = 1,
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
  b1 <- c(1, 2)
  b2 <- c(2, 3)
  b3 <- c(- 1, 4)

  mcmc <- mcmc %>%
    replace_mcmc("mod", "b1[1]", b1) %>%
    replace_mcmc("mod", "b1[2]", b2) %>%
    replace_mcmc("mod", "b1[3]", b3)
  obs <- get_extreme(
    mcmc$mod,
    time = NULL,
    greater = TRUE,
    lower = lower,
    upper = upper,
    index = NULL
  )
  exp <- tibble::tibble(doses = c(doses[2], doses[3])) %>%
    dplyr::mutate(
      extreme_responses = c(b2[1], b3[2]),
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
  exp <- tibble::tibble(doses = c(doses[3], doses[1])) %>%
    dplyr::mutate(
      extreme_responses = c(b3[1], b1[2]),
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
  exp <- tibble::tibble(doses = c(doses[3], doses[1])) %>%
    dplyr::mutate(
      extreme_responses = c(b3[1], b1[2]),
      greater = FALSE
    ) %>%
    dplyr::slice(2)
  expect_equal(obs, exp)

  # change upper bound
  obs <- get_extreme(
    mcmc$mod,
    time = NULL,
    greater = TRUE,
    lower = lower,
    upper = 4,
    index = NULL
  )
  exp <- tibble::tibble(doses = c(doses[2], doses[2])) %>%
    dplyr::mutate(
      extreme_responses = c(b2[1], b2[2]),
      greater = TRUE
    )
  expect_equal(obs, exp)

  # change lower bound
  obs <- get_extreme(
    mcmc$mod,
    time = NULL,
    greater = TRUE,
    lower = 4,
    upper = upper,
    index = NULL
  )
  exp <- tibble::tibble(doses = c(doses[3], doses[3])) %>%
    dplyr::mutate(
      extreme_responses = c(b3[1], b3[2]),
      greater = TRUE
    )
  expect_equal(obs, exp)
})

test_that("get_extreme.independent() longitudinal", {
  times <- c(0, 10)
  t_max <- max(times)
  doses <- c(0, 2, 4)
  data <- dreamer_data_linear(
    n_cohorts = c(10, 25, 30),
    dose = doses,
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
    mod = model_independent(
      mu_b1 = 0,
      sigma_b1 = 1,
      shape = 1,
      rate = .01,
      longitudinal = model_longitudinal_linear(0, 1, t_max)
    )
  )

  lower <- min(attr(mcmc, "doses"))
  upper <- max(attr(mcmc, "doses"))

  a <- c(.1, .2)
  b1 <- c(1, 2)
  b2 <- c(2, 3)
  b3 <- c(- 1, 4)

  mcmc <- mcmc %>%
    replace_mcmc("mod", "a", a) %>%
    replace_mcmc("mod", "b1[1]", b1) %>%
    replace_mcmc("mod", "b1[2]", b2) %>%
    replace_mcmc("mod", "b1[3]", b3)

  time <- 3
  obs <- get_extreme(
    mcmc$mod,
    time = 3,
    greater = TRUE,
    lower = lower,
    upper = upper,
    index = NULL
  )
  exp <- tibble::tibble(doses = c(doses[2], doses[3])) %>%
    dplyr::mutate(
      extreme_responses = a + time / t_max * c(b2[1], b3[2]),
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
  exp <- tibble::tibble(doses = c(doses[3], doses[1])) %>%
    dplyr::mutate(
      extreme_responses = a + time / t_max * c(b3[1], b1[2]),
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
  exp <- tibble::tibble(doses = c(doses[3], doses[1])) %>%
    dplyr::mutate(
      extreme_responses = a + time / t_max * c(b3[1], b1[2]),
      greater = FALSE
    ) %>%
    dplyr::slice(2)
  expect_equal(obs, exp)

  # change upper bound
  obs <- get_extreme(
    mcmc$mod,
    time = time,
    greater = TRUE,
    lower = lower,
    upper = 3,
    index = NULL
  )
  exp <- tibble::tibble(doses = c(doses[2], doses[2])) %>%
    dplyr::mutate(
      extreme_responses = a + time / t_max * c(b2[1], b2[2]),
      greater = TRUE
    )
  expect_equal(obs, exp)

  # change lower bound
  obs <- get_extreme(
    mcmc$mod,
    time = time,
    greater = TRUE,
    lower = 4,
    upper = upper,
    index = NULL
  )
  exp <- tibble::tibble(doses = c(doses[3], doses[3])) %>%
    dplyr::mutate(
      extreme_responses = a + time / t_max * c(b3[1], b3[2]),
      greater = TRUE
    )
  expect_equal(obs, exp)
})
