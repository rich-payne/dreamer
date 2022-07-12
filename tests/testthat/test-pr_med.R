test_that("pr_med() linear", {
  data <- dreamer_data_linear(
    n_cohorts = c(10, 25, 30),
    dose = c(0, 2, 4),
    b1 = .5,
    b2 = 3,
    sigma = .5
  )
  out <- dreamer_mcmc(
    data = data,
    n_iter = 2,
    n_chains = 2,
    convergence_warn = FALSE,
    silent = TRUE,
    mod = model_linear(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      shape = 1,
      rate = .01
    )
  )
  b1 <- c(1, 5,  1,  5)
  b2 <- c(1, 1, - 1, - 1)

  out_new <- out %>%
    replace_mcmc("mod", "b1", b1) %>%
    replace_mcmc("mod", "b2", b2)

  expect_equal(pr_med(out_new, csd = 0)$prob, c(1, 0, 0))
  expect_equal(pr_med(out_new, csd = 1.5)$prob, c(.5, .25, 0))
  expect_equal(
    pr_med(out_new, csd = 0, greater = FALSE)$prob,
    c(0, .25, 0)
  )
  expect_equal(
    pr_med(out_new, csd = 1.5, greater = FALSE)$prob,
    c(.5, 0, .25)
  )

  expect_equal(pr_med(out_new, csd = 0, reference_dose = 1)$prob, c(.5, .5, 0))
  expect_equal(pr_med(out_new, csd = 1.5, reference_dose = 1)$prob, c(0, 0, .5))
  expect_equal(
    pr_med(out_new, csd = 0, reference_dose = 1, greater = FALSE)$prob,
    c(.5, .5, 0)
  )
  expect_equal(
    pr_med(out_new, csd = 1.5, reference_dose = 1, greater = FALSE)$prob,
    c(1, 0, 0)
  )
})


test_that("pr_med() linear, longitudinal", {
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
  out <- dreamer_mcmc(
    data = data,
    n_iter = 2,
    n_chains = 2,
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
  a <- c(.1, .2, .1, .2)
  b1 <- c(1, 5,  1,  5)
  b2 <- c(1, 1, - 1, - 1)

  out_new <- out %>%
    replace_mcmc("mod", "a", a) %>%
    replace_mcmc("mod", "b1", b1) %>%
    replace_mcmc("mod", "b2", b2)


  expect_equal(pr_med(out_new, csd = 0, time = 4)$prob, c(1, 0, 0))
  expect_equal(
    pr_med(out_new, csd = 1.5, time = 4)$prob,
    c(.5, 0, .25)
  )
  expect_equal(
    pr_med(out_new, csd = 0, time = 4, greater = FALSE)$prob,
    c(0, .25, 0)
  )
  expect_equal(
    pr_med(out_new, csd = 1.5, time = 4, greater = FALSE)$prob,
    c(.5, .25, 0)
  )

  expect_equal(
    pr_med(out_new, csd = 0, time = 4, reference_dose = 1)$prob,
    c(.5, .5, 0)
  )
  expect_equal(
    pr_med(out_new, csd = 1.5, time = 4, reference_dose = 1)$prob,
    c(0, 0, 0)
  )
  expect_equal(
    pr_med(
      out_new,
      csd = 0,
      reference_dose = 1,
      time = 4,
      greater = FALSE
    )$prob,
    c(.5, .5, 0)
  )
  expect_equal(
    pr_med(
      out_new,
      csd = 1.5,
      reference_dose = 1,
      time = 4,
      greater = FALSE
    )$prob,
    c(1, 0, 0)
  )
})
