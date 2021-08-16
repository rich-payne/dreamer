test_that("misc helpers in dreamer.R", {
  suppressMessages(
    expect_message(
      check_aggregate(data.frame(sample_var = 1)),
      class = "dreamer"
    )
  )
})
