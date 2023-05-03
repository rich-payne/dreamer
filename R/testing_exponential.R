
library(tidyverse)
library(DoseFinding)

#exponential example
test_example <- dreamer_data_independent(
  n_cohorts = rep(50, 5),
  doses = c(0, 5, 15, 25, 50),
  b1 = c(-100, -99, -97, -95, -50),
  sigma = 50
)

independent <- dreamer_mcmc(
  data = test_example,
  independent = model_independent(
    mu_b1 = -100,
    sigma_b1 = 1000,
    shape = 1,
    rate = 0.001,
    w_prior = 1
  )
)

plot(independent, data = test_example)

exponential <- dreamer_mcmc(
  data = test_example,
  exp = model_exp(
    mu_b1 = -100,
    sigma_b1 = 1000,
    mu_b2 = 0,
    sigma_b2 = 10,
    mu_b3 = 0,
    sigma_b3 = 10,
    shape = 1,
    rate = 0.001,
    w_prior = 1
  ),
  n_burn = 1e4,
  n_iter = 1e4,
  n_chains = 3
)

plot(exponential, data = test_example)
summary(exponential)

#bma example
bma <- dreamer_mcmc(
  data = test_example,
  exp = model_exp(
    mu_b1 = -100,
    sigma_b1 = 1000,
    mu_b2 = 0,
    sigma_b2 = 10,
    mu_b3 = 0,
    sigma_b3 = 1,
    shape = 1,
    rate = 0.001,
    w_prior = 1/3
  ),
  linear = model_linear(
    mu_b1 = -100,
    sigma_b1 = 1000,
    mu_b2 = 0,
    sigma_b2 = 1000,
    shape = 1,
    rate = 0.001,
    w_prior = 1/3
  ),
  emax = model_emax(
    mu_b1 = -100,
    sigma_b1 = 1000,
    mu_b2 = 0,
    sigma_b2 = 1000,
    mu_b3 = log(15),
    sigma_b3 = log(15),
    mu_b4 = 1,
    sigma_b4 = 2,
    shape = 1,
    rate = 0.001,
    w_prior = 1/3
  ),
  
  n_burn = 1e4,
  n_iter = 1e4,
  n_chains = 3
  
)

plot(bma, data = test_example)
plot(bma$exp, data = test_example)
summary(bma)

# try starting values from dose-finding package
test_example_df <- dreamer_data_independent(
  n_cohorts = rep(50, 5),
  doses = c(0, 5, 15, 25, 50),
  b1 = c(-100, -99, -97, -95, -50),
  sigma = 50
)

df <- fitMod(dose = test_example_df$dose, resp = test_example_df$response, model = "exponential")
delta_start <- -1/df$coefs[3]

exponential <- dreamer_mcmc(
  data = test_example_df,
  exp = model_exp(
    mu_b1 = -100,
    sigma_b1 = 1000,
    mu_b2 = 0,
    sigma_b2 = 10,
    mu_b3 = delta_start,
    sigma_b3 = 1,
    shape = 1,
    rate = 0.001,
    w_prior = 1
  ),
  n_burn = 1e4,
  n_iter = 1e4,
  n_chains = 3
)

plot(exponential, data = test_example_df)
summary(exponential)
