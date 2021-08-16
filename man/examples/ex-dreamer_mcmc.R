set.seed(888)
data <- dreamer_data_linear(
  n_cohorts = c(20, 20, 20),
  dose = c(0, 3, 10),
  b1 = 1,
  b2 = 3,
  sigma = 5
)

# Bayesian model averaging
output <- dreamer_mcmc(
  data = data,
  n_adapt = 1e3,
  n_burn = 1e2,
  n_iter = 1e3,
  n_chains = 2,
  silent = TRUE,
  mod_linear = model_linear(
    mu_b1 = 0,
    sigma_b1 = 1,
    mu_b2 = 0,
    sigma_b2 = 1,
    shape = 1,
    rate = .001,
    w_prior = 1 / 2
  ),
  mod_quad = model_quad(
    mu_b1 = 0,
    sigma_b1 = 1,
    mu_b2 = 0,
    sigma_b2 = 1,
    mu_b3 = 0,
    sigma_b3 = 1,
    shape = 1,
    rate = .001,
    w_prior = 1 / 2
  )
)
# posterior weights
output$w_post
# plot posterior dose response
plot(output)

# LONGITUDINAL
library(ggplot2)
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

output_long <- dreamer_mcmc(
   data = data_long,
   n_adapt = 1e3,
   n_burn = 1e2,
   n_iter = 1e3,
   n_chains = 2,
   silent = TRUE, # make rjags be quiet,
   mod_linear = model_linear(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      shape = 1,
      rate = .001,
      w_prior = 1 / 2, # prior probability of the model
      longitudinal = model_longitudinal_itp(
         mu_a = 0,
         sigma_a = 1,
         a_c1 = 0,
         b_c1 = 1,
         t_max = 52
      )
   ),
   mod_quad = model_quad(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      mu_b3 = 0,
      sigma_b3 = 1,
      shape = 1,
      rate = .001,
      w_prior = 1 / 2,
      longitudinal = model_longitudinal_linear(
         mu_a = 0,
         sigma_a = 1,
         t_max = 52
      )
   )
)

\dontrun{
# plot longitudinal dose-response profile
plot(output_long, data = data_long)
plot(output_long$mod_quad, data = data_long) # single model

# plot dose response at final timepoint
plot(output_long, data = data_long, times = 52)
plot(output_long$mod_quad, data = data_long, times = 52) # single model
}