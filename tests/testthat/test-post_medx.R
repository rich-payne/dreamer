test_that("post_medx", {
  n_chains <- 2
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
    n_iter = 5,
    silent = TRUE,
    convergence_warn = FALSE,
    n_chains = n_chains
  )

  small_bound <- .05
  lower <- 0
  upper <- 5
  b1 <- 1:10
  b2 <- c(-c(5:1), 1:5)
  mcmc <- mcmc %>%
    replace_mcmc("mod", "b1", b1) %>%
    replace_mcmc("mod", "b2", b2)
  max_dose_gt <- c(rep(lower, 5), rep(upper, 5))
  ed100_gt <- b1 + b2 * max_dose_gt
  ed50_gt <- .5 * (ed100_gt - small_bound) + small_bound
  ed50_dose_gt <- (ed50_gt - b1) / b2
  ed50_dose_gt[
    !(ed50_dose_gt >= lower & ed50_dose_gt <= upper) |
      max_dose_gt < ed50_dose_gt
  ] <- NA

  post_medx(
    mcmc,
    ed = 50,
    greater = TRUE,
    lower = lower,
    upper = upper,
    small_bound = small_bound
  )$stats %>%
    dplyr::slice(1) %>%
    as.numeric() %>%
    expect_equal(
      c(
        50,
        mean(!is.na(ed50_dose_gt)),
        mean(ed50_dose_gt, na.rm = TRUE),
        quantile(
          ed50_dose_gt,
          prob = c(.025, .975),
          na.rm = TRUE,
          names = FALSE
        )
      )
   )

  max_dose_lt <- c(rep(upper, 5), rep(lower, 5))
  ed100_lt <- b1 + b2 * max_dose_lt
  ed50_lt <- .5 * (ed100_lt - small_bound) + small_bound
  ed50_dose_lt <- (ed50_lt - b1) / b2
  ed50_dose_lt[
    !(ed50_dose_lt >= lower & ed50_dose_lt <= upper) |
      max_dose_lt < ed50_dose_lt
  ] <- NA

  post_medx(
    mcmc,
    ed = 50,
    greater = FALSE,
    lower = lower,
    upper = upper,
    small_bound = small_bound
  )$stats %>%
    dplyr::slice(1) %>%
    as.numeric() %>%
    expect_equal(
      c(
        50,
        mean(!is.na(ed50_dose_lt)),
        mean(ed50_dose_lt, na.rm = TRUE),
        quantile(
          ed50_dose_lt,
          prob = c(.025, .975),
          na.rm = TRUE,
          names = FALSE
        )
      )
    )
})

test_that("post_medx longitudinal", {
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
      longitudinal = model_longitudinal_linear(0, 1, t_max)
    ),
    n_iter = 10,
    silent = TRUE,
    convergence_warn = FALSE,
    n_chains = 1
  )

  small_bound <- .05
  lower <- 0
  upper <- 5
  a <- c(0:9) / 10
  b1 <- 1:10
  b2 <- c(-c(5:1), 1:5)
  mcmc <- mcmc %>%
    replace_mcmc("mod", "a", a) %>%
    replace_mcmc("mod", "b1", b1) %>%
    replace_mcmc("mod", "b2", b2)
  time <- 3
  max_dose_gt <- c(rep(lower, 5), rep(upper, 5))
  ed100_gt <- a + time / t_max * (b1 + b2 * max_dose_gt)
  ed50_gt <- .5 * (ed100_gt - small_bound) + small_bound
  ed50_dose_gt <- (ed50_gt - a - b1 * time / t_max) / (b2 * time / t_max)
  ed50_dose_gt[
    !(ed50_dose_gt >= lower & ed50_dose_gt <= upper) |
      max_dose_gt < ed50_dose_gt
  ] <- NA

  post_medx(
    mcmc,
    ed = 50,
    time = time,
    greater = TRUE,
    lower = lower,
    upper = upper,
    small_bound = small_bound
  )$stats %>%
    dplyr::slice(1) %>%
    as.numeric() %>%
    expect_equal(
      c(
        50,
        mean(!is.na(ed50_dose_gt)),
        mean(ed50_dose_gt, na.rm = TRUE),
        quantile(
          ed50_dose_gt,
          prob = c(.025, .975),
          na.rm = TRUE,
          names = FALSE
        )
      )
    )

  max_dose_lt <- c(rep(upper, 5), rep(lower, 5))
  ed100_lt <- a + time / t_max * (b1 + b2 * max_dose_lt)
  ed50_lt <- .5 * (ed100_lt - small_bound) + small_bound
  ed50_dose_lt <- (ed50_lt - a - b1 * time / t_max) / (b2 * time / t_max)
  ed50_dose_lt[
    !(ed50_dose_lt >= lower & ed50_dose_lt <= upper) |
      max_dose_lt < ed50_dose_lt
  ] <- NA

  post_medx(
    mcmc,
    ed = 50,
    time = time,
    greater = FALSE,
    lower = lower,
    upper = upper,
    small_bound = small_bound
  )$stats %>%
    dplyr::slice(1) %>%
    as.numeric() %>%
    expect_equal(
      c(
        50,
        mean(!is.na(ed50_dose_lt)),
        mean(ed50_dose_lt, na.rm = TRUE),
        quantile(
          ed50_dose_lt,
          prob = c(.025, .975),
          na.rm = TRUE,
          names = FALSE
        )
      )
    )
})

test_that("post_medx.dreamer_bma()", {
  data <- dreamer_data_linear(n_cohorts = c(10, 20, 30), c(1, 3, 5), 1, 2, 2)
  set.seed(88332)
  mcmc <- dreamer_mcmc(
    data,
    mod = model_linear(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      shape = 1,
      rate = .01,
      w_prior = .5
    ),
    quad = model_quad(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      mu_b3 = 0,
      sigma_b3 = 1,
      shape = 1,
      rate = .01,
      w_prior = .5
    ),
    n_iter = 5,
    silent = TRUE,
    convergence_warn = FALSE,
    n_chains = 1
  )

  model_index <- attr(mcmc, "model_index")
  exp <- dplyr::bind_rows(
    post_medx(mcmc$mod, ed = 50, return_samples = TRUE)$samps %>%
      dplyr::slice(which(model_index == 1)),
    post_medx(mcmc$quad, ed = 50, return_samples = TRUE)$samps %>%
      dplyr::slice(which(model_index == 2))
  ) %>%
    dplyr::group_by(ed) %>%
    dplyr::summarize(
      pr_edx_exists = mean(!is.na(dose)),
      mean = mean(dose),
      `2.5%` = quantile(dose, prob = .025),
      `97.5%` = quantile(dose, prob = .975),
      .groups = "drop"
    ) %>%
    as.data.frame()
  expect_equal(post_medx(mcmc, ed = 50)$stats, exp)
})

test_that("post_medx() errors for independent models", {
  data <- dreamer_data_linear(n_cohorts = c(10, 20, 30), c(1, 3, 5), 1, 2, 2)
  mcmc <- dreamer_mcmc(
    data,
    mod = model_independent(0, 1, 1, 1),
    n_iter = 5,
    silent = TRUE,
    convergence_warn = FALSE,
    n_chains = 1
  )
  expect_error(post_medx(mcmc), class = "dreamer")

  data <- dreamer_data_linear_binary(
    n_cohorts = c(10, 20, 30),
    dose = c(1, 3, 5),
    b1 = 1,
    b2 = 2,
    link = "logit"
  )
  mcmc <- dreamer_mcmc(
    data,
    mod = model_independent_binary(0, 1, link = "logit"),
    n_iter = 5,
    silent = TRUE,
    convergence_warn = FALSE,
    n_chains = 1
  )
  expect_error(post_medx(mcmc), class = "dreamer")
})
