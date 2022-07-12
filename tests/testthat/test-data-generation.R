# setup
test_data <- function(func, longitudinal, n_rows, ...) {
  d <- rlang::exec(
    func,
    longitudinal = longitudinal,
    ...
  )
  testthat_data(d, longitudinal, n_rows)
}

testthat_data <- function(dat, longitudinal, n_rows) {
  expect_true(is.data.frame(dat))
  expect_true(
    all(rlang::has_name(dat, c("cohort", "dose", "subject", "response")))
  )
  if (!is.null(longitudinal)) {
    expect_true(rlang::has_name(dat, c("time")))
  }
  expect_equal(nrow(dat), n_rows)
}

# longitudinal models
long_models <- c(list(NULL), "linear", "itp", "idp")
t_max <- 5
long_params <- tibble::tibble(
  b1 = .1,
  b2 = .5,
  b3 = .7,
  b4 = 1,
  sigma = 1.2,
  a = 2,
  c1 = .1,
  c2 = - 1,
  gam = .75,
  d1 = 1,
  d2 = 4,
  times = list(c(1, 3, 5)),
  doses = list(c(1:3)),
  scale = 1.2 * max(doses[[1]]),
  n_cohorts = list(c(10, 15, 10)),
  t_max = !!t_max
)

test_that("data generation: continuous", {
  dat <- tidyr::expand_grid(
    dose_model = c(
      "linear", "quad", "loglinear", "logquad", "emax", "exp", "beta"
    ),
    longitudinal = !!long_models,
    long_params
  ) %>%
    dplyr::mutate(
      func = paste0("dreamer_data_", dose_model),
      n_rows = dplyr::case_when(
        vapply(longitudinal, is.null, logical(1)) ~ 35,
        TRUE ~ 3 * 35
      )
    )
  purrr::pmap(
    dplyr::select(dat, - dose_model),
    test_data
  )
})

test_that("data generation: continuous independent", {
  dat_ind <- tidyr::expand_grid(
    dose_model = "independent",
    longitudinal = !!long_models,
    long_params
  ) %>%
    dplyr::mutate(
      func = paste0("dreamer_data_", dose_model),
      b1 = list(c(.5, 1.5, 2.5)),
      n_rows = dplyr::case_when(
        vapply(longitudinal, is.null, logical(1)) ~ 35,
        TRUE ~ 3 * 35
      )
    )
  purrr::pmap(
    dplyr::select(dat_ind, - dose_model),
    test_data
  )
})

test_that("data generation: binary", {
  dat_bin <- tidyr::expand_grid(
    dose_model = c(
      "linear", "quad", "loglinear", "logquad", "emax", "exp", "beta"
    ),
    link = c("logit", "probit"),
    longitudinal = !!long_models,
    long_params
  ) %>%
    dplyr::mutate(
      func = paste0("dreamer_data_", dose_model, "_binary"),
      n_rows = dplyr::case_when(
        vapply(longitudinal, is.null, logical(1)) ~ 35,
        TRUE ~ 3 * 35
      )
    )
  purrr::pmap(
    dplyr::select(dat_bin, - dose_model, - sigma),
    test_data
  )
})

test_that("data generation: binary independent", {
  dat_bin_ind <- tidyr::expand_grid(
    dose_model = "independent_binary",
    longitudinal = !!long_models,
    link = c("logit", "probit"),
    long_params
  ) %>%
    dplyr::mutate(
      func = paste0("dreamer_data_", dose_model),
      b1 = list(c(.5, 1.5, 2.5)),
      n_rows = dplyr::case_when(
        vapply(longitudinal, is.null, logical(1)) ~ 35,
        TRUE ~ 3 * 35
      )
    )
  purrr::pmap(
    dplyr::select(dat_bin_ind, - dose_model, - sigma),
    test_data
  )
})
