test_that("misc helper functions", {
  expect_equal(solve_quad(1, .01, 1), rep(NA_real_, 2))
  expect_equal(expand_mean_time(1:2, matrix(1, 1, 1)), 1:2)
  expect_equal(get_min_dose(matrix(NA, 2, 5)), rep(NA_real_, 5))
  expect_equal(get_min_na(numeric(0)), NA_real_)
})
