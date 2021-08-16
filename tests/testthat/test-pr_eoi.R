rep_true <- function(x, doses, dose) {
  out <- rep(x, length(doses))
  out[which(dose == doses)] <- 0
  return(out)
}

check_pr_eoi <- function(mcmc, reference_dose = NULL) {
  # reference_dose is NULL
  post <- posterior(mcmc, probs = c(.25, .75))$stats
  n_mcmc <- nrow(as.matrix(mcmc[[4]]))
  tol <- 1 / n_mcmc + .Machine$double.eps
  pr_eoi(mcmc, eoi = post$`25.00%`, dose = post$dose)$prob %>%
    expect_equal(rep(.75, nrow(post)), tolerance = tol)
  pr_eoi(mcmc, eoi = post$`75.00%`, dose = post$dose)$prob %>%
    expect_equal(rep(.25, nrow(post)), tolerance = tol)

  # reference_dose
  post2 <- posterior(
    mcmc,
    probs = c(.25, .75),
    reference_dose = reference_dose
  )$stats
  pr_eoi(
    mcmc,
    eoi = post$`25.00%`,
    dose = post$dose,
    reference_dose = reference_dose
  )$prob %>%
    expect_equal(rep_true(.75, post$dose, reference_dose), tolerance = tol)
  pr_eoi(
    mcmc,
    eoi = post$`75.00%`,
    dose = post$dose,
    reference_dose = reference_dose
  )$prob %>%
    expect_equal(rep_true(.25, post$dose, reference_dose), tolerance = tol)
}

check_pr_eoi_long <- function(mcmc, times = c(1, 5), reference_dose = NULL) {
  # assume arm 2 is NULL
  post <- posterior(mcmc, probs = c(.25, .75), times = times)$stats
  n_mcmc <- nrow(as.matrix(mcmc[[4]]))
  tol <- 1 / n_mcmc + .Machine$double.eps
  purrr::map(
    times,
    function(time) {
      xx <- dplyr::filter(post, time == !!time)
      pr_eoi(mcmc, eoi = xx$`25.00%`, dose = xx$dose, time = time)$prob %>%
        expect_equal(rep(.75, nrow(xx)), tolerance = tol)
      pr_eoi(mcmc, eoi = xx$`75.00%`, dose = xx$dose, time = time)$prob %>%
        expect_equal(rep(.25, nrow(xx)), tolerance = tol)
    }
  )
  # reference_dose
  post2 <- posterior(
    mcmc,
    probs = c(.25, .75),
    reference_dose = reference_dose,
    times = times
  )$stats
  purrr::map(
    times,
    function(time) {
      xx <- dplyr::filter(post2, time == !!time)
      pr_eoi(
        mcmc,
        eoi = xx$`25.00%`,
        dose = xx$dose,
        reference_dose = reference_dose,
        time = time
      )$prob %>%
        expect_equal(rep_true(.75, xx$dose, reference_dose), tolerance = tol)
      pr_eoi(
        mcmc,
        eoi = xx$`75.00%`,
        dose = xx$dose,
        reference_dose = reference_dose,
        time = time
      )$prob %>%
        expect_equal(rep_true(.25, xx$dose, reference_dose), tolerance = tol)
    }
  )
}

test_that("pr_eoi()", {
  data <- dreamer_data_linear(
    n_cohorts = c(10, 25, 30),
    dose = c(0, 2, 4),
    b1 = .5,
    b2 = 3,
    sigma = .5
  )
  out <- dreamer_mcmc(
    data = data,
    n_iter = 1000,
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
  check_pr_eoi(out, reference_dose = 1.5)
})

test_that("pr_eoi.dreamer() longitudinal", {
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
  check_pr_eoi_long(out, times = c(0, 5, 10), reference_dose = 1.5)
})

test_that("pr_eoi() with grid", {
  data <- dreamer_data_linear(
    n_cohorts = c(10, 25, 30),
    dose = c(0, 2, 4),
    b1 = .5,
    b2 = .3,
    sigma = .5
  )
  out <- dreamer_mcmc(
    data = data,
    n_iter = 1000,
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
  obs <- pr_eoi(out, eoi = .75, dose = c(1, 1.1, 1.2))
  obs2 <- pr_eoi(out, eoi = rep(.75, 3), dose = c(1, 1.1, 1.2))
  obs3 <- pr_eoi(out, eoi = c(.75, .76), dose = c(1, 1.1, 1.2))
  exp <- dplyr::bind_rows(
    pr_eoi(out, eoi = .75, dose = 1),
    pr_eoi(out, eoi = .75, dose = 1.1),
    pr_eoi(out, eoi = .75, dose = 1.2)
  )
  exp3 <- dplyr::bind_rows(
    exp,
    pr_eoi(out, eoi = .76, dose = 1),
    pr_eoi(out, eoi = .76, dose = 1.1),
    pr_eoi(out, eoi = .76, dose = 1.2)
  )
  expect_equal(obs, exp)
  expect_equal(obs2, exp)
  expect_equal(obs3, exp3)
})
