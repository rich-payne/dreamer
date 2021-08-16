test_that("pextreme()", {
  responses <- data.frame(
    dose = 0:1,
    mean_response = c(1, 2, 3, 2),
    iter = rep(1:2, each = 2)
  )

  obs <- pextreme(responses, greater = TRUE)
  exp <- tibble::tibble(
    doses = c(1, 0),
    extreme_responses = 2:3,
    greater = TRUE
  )
  expect_equal(obs, exp)

  obs <- pextreme(responses, greater = FALSE)
  exp <- tibble::tibble(
    doses = c(0, 1),
    extreme_responses = 1:2,
    greater = FALSE
  )
  expect_equal(obs, exp)
})

test_that("get_extreme.R helpers work", {
  expect_error(get_extreme(1), class = "dreamer")
})
