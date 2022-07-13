set.seed(888)
data <- dreamer_data_linear(
  n_cohorts = c(20, 20, 20),
  dose = c(0, 3, 10),
  b1 = 1,
  b2 = 3,
  sigma = 5
)

head(data)

plot(data$dose, data$response)
abline(a = 1, b = 3)

# longitudinal data
set.seed(889)
data_long <- dreamer_data_linear(
  n_cohorts = c(10, 10, 10, 10), # number of subjects in each cohort
  doses = c(.25, .5, .75, 1.5), # dose administered to each cohort
  b1 = 0, # intercept
  b2 = 2, # slope
  sigma = .5, # standard deviation,
  longitudinal = "itp",
  times = c(0, 12, 24, 52),
  t_max = 52, # maximum time
  a = .5,
  c1 = .1
)

\dontrun{
  ggplot(data_long, aes(time, response, group = dose, color = factor(dose))) +
    geom_point()
}
