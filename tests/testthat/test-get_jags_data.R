test_that("get_jags_data hyperparameters for independent model", {
  doses <- c(1, 5)
  data <- dreamer_data_linear(n_cohorts = c(10, 30), doses, 1, 2, 2)
  mod <- model_independent(
    mu_b1 = c(1, 3),
    sigma_b1 = c(.5, .75),
    shape = 1,
    rate = .01,
    doses = c(5, 1)
  )
  jags_data <- get_jags_data(mod, data, doses = doses)
  expect_equal(jags_data$mu_b1, c(3, 1))
  expect_equal(jags_data$sigma_b1, c(.75, .5))
})

test_that("get_jags_data helpers work", {
  jags_data <- data.frame(a = 1)
  expect_equal(set_scale(jags_data, data = NULL), jags_data)
  expect_null(prep_binary_jags_data(data = NULL))
})

test_that("prep_jags returns list whn data is NULL", {
  mod <- model_independent(
    mu_b1 = c(1, 3),
    sigma_b1 = c(.5, .75),
    shape = 1,
    rate = .01,
    doses = c(5, 1)
  )
  jags_data <- get_jags_data(mod, data = NULL)
  expect_true(is.list(jags_data))
  expect_equal(jags_data$n_doses, 2)

  mod <- model_independent_binary(
    mu_b1 = c(1, 3),
    sigma_b1 = c(.5, .75),
    link = "logit",
    doses = c(5, 1)
  )
  jags_data <- get_jags_data(mod, data = NULL)
  expect_true(is.list(jags_data))
  expect_equal(jags_data$n_doses, 2)
})

test_that("prep_cont_jags_data runs for aggregated data", {
  data <- data.frame(
    dose = 1:3,
    response = 11:13,
    sample_var = c(1.1, 2.2, 3.3),
    n = c(10, 20, 30)
  )
  mod <- model_independent(
    mu_b1 = c(1, 3),
    sigma_b1 = c(.5, .75),
    shape = 1,
    rate = .01,
    doses = c(5, 1)
  )
  expect_true(is.list(prep_cont_jags_data(mod, data)))
})

test_that("remove_t_max is the identity function for IDP", {
  mod <- model_independent(
    mu_b1 = c(1, 3),
    sigma_b1 = c(.5, .75),
    shape = 1,
    rate = .01,
    doses = c(5, 1),
    longitudinal = model_longitudinal_idp(1, 1, t_max = 5)
  )
  expect_equal(remove_t_max(mod$longitudinal, mod$longitudinal)$t_max, 5)
})

test_that("aggregate_binary returns data when data is NULL", {
  expect_null(aggregate_binary(NULL))
})

test_that("get_time_var()", {
  expect_equal(get_time_var(list(longitudinal = "itp")), "time")
})
