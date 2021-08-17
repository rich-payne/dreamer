set.seed(8838383)
data <- dreamer_data_linear(n_cohorts = c(10, 20, 30), c(1, 3, 5), 1, 2, 2)
mcmc <- dreamer_mcmc(
  data,
  lin = model_linear(
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
  n_iter = 100,
  n_chains = 1,
  silent = TRUE,
  convergence_warn = FALSE
)

times <- c(0, 10)
t_max <- max(times)
data_long <- dreamer_data_linear(
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

mcmc_long <- dreamer_mcmc(
  data_long,
  lin = model_linear(
    mu_b1 = 0,
    sigma_b1 = 1,
    mu_b2 = 0,
    sigma_b2 = 1,
    shape = 1,
    rate = .01,
    w_prior = .5,
    longitudinal = model_longitudinal_linear(0, 1, t_max = t_max)
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
    w_prior = .5,
    longitudinal = model_longitudinal_linear(0, 1, t_max = t_max)
  ),
  n_iter = 100,
  n_chains = 1,
  silent = TRUE,
  convergence_warn = FALSE
)

gg_save <- function(plot, filename, ...) {
  fs::dir_create("figs")
  filename <- fs::path("figs", filename)
  ggplot2::ggsave(filename, plot, width = 7, height = 7, ...)
  filename
}

test_that("plot works", {
  skip_on_cran()
  #skip_on_ci()
  out <- plot(mcmc, n_smooth = 5)
  expect_snapshot_file(gg_save(out, "plot_works.png"))
})

test_that("plot.dreamer works", {
  skip_on_cran()
  #skip_on_ci()
  out <- plot(mcmc$lin, n_smooth = 5)
  expect_snapshot_file(gg_save(out, "plot_dreamer_works.png"))
})

test_that("plot with data", {
  skip_on_cran()
  #skip_on_ci()
  out <- plot(mcmc, data = data, n_smooth = 5)
  expect_snapshot_file(gg_save(out, "plot_with_data.png"))
  # with aggregated data
  data_sum <- data %>%
    dplyr::group_by(dose) %>%
    dplyr::summarize(response = mean(response), n = n(), .groups = "drop")
  out <- plot(mcmc, data = data_sum, n_smooth = 5)
  expect_snapshot_file(gg_save(out, "plot_with_data_sum.png"))
})

test_that("predictive plots", {
  skip_on_cran()
  #skip_on_ci()
  out <- plot(mcmc, predictive = 10, n_smooth = 5)
  expect_snapshot_file(gg_save(out, "predictive_plot.png"))
})

test_that("longitudinal plots", {
  skip_on_cran()
  #skip_on_ci()
  out <- plot(mcmc_long$lin, times = c(0, t_max), n_smooth = 5)
  expect_snapshot_file(gg_save(out, "longitudinal.png"))
  out <- plot(
    mcmc_long,
    times = c(0, t_max),
    predictive = 10,
    data = data_long,
    n_smooth = 5
  )
  expect_snapshot_file(gg_save(out, "longitudinal_data_predictive.png"))
})

test_that("plots comparison", {
  skip_on_cran()
  #skip_on_ci()
  out <- plot_comparison(mcmc, n_smooth = 5)
  expect_snapshot_file(gg_save(out, "compare.png"))
  out <- plot_comparison(mcmc, data = data, n_smooth = 5)
  expect_snapshot_file(gg_save(out, "compare_data.png"))
  out <- plot_comparison(mod1 = mcmc$lin, mod2 = mcmc$quad, n_smooth = 5)
  expect_snapshot_file(gg_save(out, "compare_custom.png"))
  out <- plot_comparison(mcmc_long, n_smooth = 5)
  expect_snapshot_file(gg_save(out, "compare_longitudinal.png"))
  out <- plot_comparison(mcmc_long, times = 1:3, doses = 2, n_smooth = 5)
  expect_snapshot_file(gg_save(out, "compare_longitudinal_mulitiple_times.png"))
  out <- plot_comparison(mod1 = mcmc_long$lin, doses = 2, n_smooth = 5)
  expect_snapshot_file(
    suppressMessages(gg_save(out, "compare_longitudinal_single_time.png"))
  )
})

test_that("traceplots work", {
  skip_on_cran()
  #skip_on_ci()
  png("figs/traceplots.png")
    plot_trace(mcmc)
  dev.off()
  expect_snapshot_file("figs/traceplots.png")
  png("figs/traceplots_single.png")
    plot_trace(mcmc$lin)
  dev.off()
  expect_snapshot_file("figs/traceplots_single.png")
})

test_that("independent predictive plot", {
  skip_on_cran()
  #skip_on_ci()
  data <- dreamer_data_linear(n_cohorts = c(10, 20, 30), c(1, 3, 5), 1, 2, 2)
  mcmc <- dreamer_mcmc(
    data,
    lin = model_independent(0, 10, 1, .01),
    n_iter = 100,
    n_chains = 1,
    silent = TRUE,
    convergence_warn = FALSE
  )
  out <- plot(mcmc, predictive = 10, data = data, n_smooth = 5)
  expect_snapshot_file(gg_save(out, "independent_predictive.png"))
})

test_that("dose adjusted plot", {
  skip_on_cran()
  #skip_on_ci()
  out <- plot(mcmc, reference_dose = 0, n_smooth = 5)
  expect_snapshot_file(gg_save(out, "reference_dose.png"))
})

test_that("dreamer_plot_prior", {
  skip_on_cran()
  #skip_on_ci()
  out <- dreamer_plot_prior(
    n_samples = 1e2,
    doses = 1:3,
    mod = model_linear(0, 1, 0, 1, 1, .01)
  )
  expect_snapshot_file(gg_save(out, "prior.png"))
})

test_that("dreamer_plot_prior longitudinal single timepoint", {
  skip_on_cran()
  #skip_on_ci()
  out <- dreamer_plot_prior(
    n_samples = 1e2,
    doses = 1:3,
    mod = model_linear(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      shape = 1,
      rate = .01,
      longitudinal = model_longitudinal_linear(0, 1, t_max = 5)
    ),
    times = 1
  )
  expect_snapshot_file(gg_save(out, "prior_longitudinal_single.png"))
})

test_that("dreamer_plot_prior longitudinal multiple timepoints", {
  skip_on_cran()
  #skip_on_ci()
  out <- dreamer_plot_prior(
    n_samples = 1e2,
    doses = 1:3,
    mod = model_linear(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      shape = 1,
      rate = .01,
      longitudinal = model_longitudinal_linear(0, 1, t_max = 5)
    ),
    times = 1:3
  )
  expect_snapshot_file(gg_save(out, "prior_longitudinal_multiple.png"))
})

test_that("dreamer_plot_prior plot draws", {
  skip_on_cran()
  #skip_on_ci()
  out <- dreamer_plot_prior(
    n_samples = 1e2,
    doses = 1:3,
    mod = model_linear(0, 1, 0, 1, 1, .01),
    plot_draws = TRUE
  )
  expect_snapshot_file(gg_save(out, "prior_draws.png"))
})

test_that("dreamer_plot_prior plot draws longitudinal one timepoint", {
  skip_on_cran()
  #skip_on_ci()
  out <- dreamer_plot_prior(
    n_samples = 1e2,
    doses = 1:3,
    mod = model_linear(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      shape = 1,
      rate = .01,
      longitudinal = model_longitudinal_linear(0, 1, t_max = 5)
    ),
    times = 1,
    plot_draws = TRUE
  )
  expect_snapshot_file(gg_save(out, "prior_draws_longitudinal_single.png"))
})

test_that("dreamer_plot_prior plot draws longitudinal multiple timepoints", {
  skip_on_cran()
  #skip_on_ci()
  out <- dreamer_plot_prior(
    n_samples = 1e2,
    doses = 1:3,
    mod = model_linear(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      shape = 1,
      rate = .01,
      longitudinal = model_longitudinal_linear(0, 1, t_max = 5)
    ),
    times = 1:3,
    plot_draws = TRUE
  )
  expect_snapshot_file(gg_save(out, "prior_draws_longitudinal_multiple.png"))
})

test_that("bar_width for single dose", {
  expect_equal(bar_width(1), 0.1)
})

test_that("binary plot with data", {
  skip_on_cran()
  #skip_on_ci()
  data <- dreamer_data_linear_binary(
    n_cohorts = c(10, 20, 30),
    dose = c(1, 3, 5),
    b1 = -1,
    b2 = .5,
    link = "logit"
  )
  mcmc <- dreamer_mcmc(
    data,
    lin = model_linear_binary(
      mu_b1 = 0,
      sigma_b1 = 3,
      mu_b2 = 0,
      sigma_b2 = 3,
      link = "logit"
    ),
    n_iter = 100,
    n_chains = 1,
    silent = TRUE,
    convergence_warn = FALSE
  )
  out <- plot(mcmc, data = data, n_smooth = 5)
  expect_snapshot_file(gg_save(out, "binary_data.png"))
  # with aggregated data
  data_sum <- data %>%
    dplyr::group_by(dose) %>%
    dplyr::summarize(response = sum(response), n = n(), .groups = "drop")
  out <- plot(mcmc, data = data_sum, n_smooth = 5)
  expect_snapshot_file(gg_save(out, "binary_data_sum.png"))
})
