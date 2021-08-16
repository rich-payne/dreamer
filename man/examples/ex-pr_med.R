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
 n_burn = 1e3,
 n_iter = 1e4,
 n_chains = 2,
 silent = FALSE,
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

pr_med(output, csd = 10)

# difference of two doses
pr_med(output, csd = 3, reference_dose = 0)

# single model
pr_med(output$mod_quad, csd = 10)
