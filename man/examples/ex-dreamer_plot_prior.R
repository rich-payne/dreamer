# Plot prior for one model
set.seed(8111)
dreamer_plot_prior(
 doses = c(0, 2.5, 5),
 mod_quad_binary = model_quad_binary(
   mu_b1 = -.5,
   sigma_b1 = .2,
   mu_b2 = -.5,
   sigma_b2 = .2,
   mu_b3 = .5,
   sigma_b3 = .1,
   link = "logit",
   w_prior = 1
 )
)

# plot individual draws
dreamer_plot_prior(
 doses = seq(from = 0, to = 5, length.out = 50),
 n_samples = 100,
 plot_draws = TRUE,
 mod_quad_binary = model_quad_binary(
   mu_b1 = -.5,
   sigma_b1 = .2,
   mu_b2 = -.5,
   sigma_b2 = .2,
   mu_b3 = .5,
   sigma_b3 = .1,
   link = "logit",
   w_prior = 1
 )
)

# plot prior from mixture of models
dreamer_plot_prior(
 doses = c(0, 2.5, 5),
 mod_linear_binary = model_linear_binary(
   mu_b1 = -1,
   sigma_b1 = .1,
   mu_b2 = 1,
   sigma_b2 = .1,
   link = "logit",
   w_prior = .75
 ),
 mod_quad_binary = model_quad_binary(
   mu_b1 = -.5,
   sigma_b1 = .2,
   mu_b2 = -.5,
   sigma_b2 = .2,
   mu_b3 = .5,
   sigma_b3 = .1,
   link = "logit",
   w_prior = .25
 )
)
