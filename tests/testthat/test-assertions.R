test_that("assertions", {
  expect_error(check_link("something"))
  expect_error(check_time(1:5))
  expect_error(check_ind_model_parms(doses = NULL, mu_b1 = 1, sigma_b1 = 1:2))
  expect_error(check_ind_model_parms(doses = NULL, mu_b1 = 1:2, sigma_b1 = 1:2))
  expect_error(check_ind_model_parms(doses = 1, mu_b1 = 1:2, sigma_b1 = 1:2))
  expect_error(check_ind_model_parms(doses = 1:2, mu_b1 = 1, sigma_b1 = 1))
  expect_error(check_dose_length(1:2))
  x <- 3
  attr(x, "doses") <- 1:3
  expect_error(check_dose_index(dose_index = numeric(0), dose = 3, x))
  expect_error(
    dreamer_mean_longitudinal_mcmc(longitudinal_model = "dummy_model")
  )
  expect_error(check_len(1, 1:2, "error"))
  expect_error(check_duplicate_dose(c(1, 1)))
  expect_equal(get_progress_bar(TRUE), "none")
  expect_equal(get_progress_bar(FALSE), "text")
  expect_true(
    is.list(
      get_jags_seed(jags_seed = 1:2, jags_rng = c("a", "b"), n_chains = 2)
    )
  )
  expect_error(expand_jags_rng(jags_rng = 1:2, n_chains = 3))
  expect_error(check_seed_len(jags_seed = 1, n_chains = 2))
  expect_error(get_scale(list(a = 1), data = NULL))
  expect_error(posterior("a"))
  expect_invisible(check_independent_doses(NULL, 1:5))
  expect_error(check_independent_doses(1:3, 1:2))
  expect_error(check_bounds(dose = 1, lower = 2, upper = 3))
  expect_error(check_bounds(dose = 4, lower = 2, upper = 3))
  expect_error(check_small_bound2(small_bound = 1:2, index = 1:3))
  expect_error(
    check_doses_data(doses = 1:2, data = data.frame(a = 1), u_doses_data = 2:3)
  )
  expect_error(check_colnames(data.frame(n = 1)))
  expect_error(check_independent_model(NULL, NULL))
  x2 <- 3
  attr(x2, "longitudinal_model") <- "dummy_longitudinal"
  expect_error(get_n_params_longitudinal(x2))
  expect_error(check_eoi_lengths(1, 1:2))
  expect_error(check_eoi_lengths(1:2, 1:2, 1))
  x3 <- list(times = 1:3)
  attr(x3, "times") <- 1:5
  expect_equal(get_time(x3, time = NULL), 5)
  attr(x3, "times") <- NULL
  expect_equal(get_time(x3, time = NULL), 3)
  expect_error(get_time(x3, time = 1:2))
  expect_error(check_small_bound(1:5, 4))
  suppressMessages(expect_message(mcmc_start_msg("my_mod", Sys.time(), FALSE)))
  suppressMessages(expect_message(mcmc_end_msg(Sys.time(), Sys.time(), FALSE)))
  expect_error(assert_w_prior(c(1, 1)))
  expect_error(assert_dreamer_dots(list(a = 1)))
  expect_error(
    assert_independent_dots(
      list(
        model_independent(1, 1, 1, 1),
        model_linear(1, 1, 1, 1, 1, 1)
      )
    )
  )
  expect_error(
    assert_binary_dots(
      list(
        model_linear(1, 1, 1, 1, 1, 1),
        model_linear_binary(1, 1, 1, 1, "logit")
      )
    )
  )
  expect_error(check_data(data.frame(dose = 1)))
  expect_error(check_data(data.frame(response = 1)))
  expect_warning(
    throw_convergence_warn(
      data.frame(model = "beta", param = "b1", gelman_upper = 1.5332)
    )
  )
  expect_error(name_models(list(3, model_1_numeric = 2)))
  expect_error(assert_names(list(a = 1), list(b = 1), "names mistmatch"))
  expect_error(check_binary(data.frame(response = 3)))
  expect_error(check_binary(data.frame(response = - 1)))
  expect_error(
    check_longitudinal(list(list(longitudinal = "a"), list(b = NULL)))
  )
  expect_error(
    check_longitudinal(list(list(longitudinal = "a")), data.frame(a = 1))
  )
  expect_error(check_inputs(1, 1, 1))
  expect_error(check_times(1:5, FALSE))
  expect_error(check_times(NULL, TRUE))
  expect_error(check_names(1))
  expect_error(assert_dose_len(1:2))
  expect_error(check_plot(1:2, 1))
  expect_error(check_plot(NULL, 1:2))
  expect_error(check_no_dots("fun_name", a = 3))
  expect_error(assert_data_reference_dose(1, 1), class = "dreamer")
  expect_error(assert_doses(NULL), class = "dreamer")
  expect_error(assert_no_dots("foo", foo2 = 5), class = "dreamer")
  expect_error(assert_reference_dose(1:2), class = "dreamer")
})
