# useful function to help with determining truth for tests
get_edx_true_lin <- function(
  ed,
  a = 0,
  b1,
  b2,
  doses,
  small_bound,
  greater,
  time_perc = 1
) {
  responses <- a + (b1 + b2 * doses) * time_perc
  print(responses)
  if (greater) {
    ed100 <- max(responses)
    edx <- small_bound + ed / 100 * (ed100 - small_bound)
    edx <- min(doses[responses > edx])
  } else {
    ed100 <- min(responses)
    edx <- small_bound + ed / 100 * (ed100 - small_bound)
    edx <- min(doses[responses < edx])
  }
  return(edx)
}

test_that("pr_medx.dreamer()", {
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

  expect_equal(pr_medx(out_new, ed = 50)$prob, c(.75, .25, 0))
  expect_equal(
    pr_medx(out_new, ed = 50, small_bound = 1)$prob,
    c(.25, .25, .25)
  )
  expect_equal(pr_medx(out_new, ed = 50, greater = FALSE)$prob, c(0, 0, .25))
  expect_equal(
    pr_medx(out_new, ed = 50, small_bound = 1, greater = FALSE)$prob,
    c(0, 0, .25)
  )
})

test_that("pr_medx.dreamer() longitudinal", {
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

  expect_equal(pr_medx(out_new, ed = 50, time = .2)$prob, c(1, 0, 0))
  expect_equal(
    pr_medx(out_new, ed = 50, small_bound = .1, time = .2)$prob,
    c(.75, .25, 0)
  )
  expect_equal(
    pr_medx(out_new, ed = 50, time = .2, greater = FALSE)$prob,
    c(0, 0, 0)
  )
  expect_equal(
    pr_medx(out_new, ed = 50, small_bound = .1, greater = FALSE)$prob,
    c(0, 0, .25)
  )
})
