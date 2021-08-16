test_that("Draw values from prior", {
  out <- dreamer_mcmc(
    data = NULL,
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
  expect_true(inherits(out, "dreamer_bma"))
})
