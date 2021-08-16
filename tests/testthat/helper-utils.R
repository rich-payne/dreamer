replace_mcmc <- function(x, mod, var, val) {
  if (is.null(val)) return(x)
  val <- assert_val(x[[mod]], val)
  n_mcmc <- length(val)
  n_chains <- length(x[[mod]])
  n_per_chain <- n_mcmc / n_chains
  for (i in seq_along(x[[mod]])) {
    index <- ((i - 1) * n_per_chain + 1):(i * n_per_chain)
    x[[mod]][[i]][, var] <- val[index]
  }
  x
}

assert_val <- function(x, val) {
  n_mcmc <- nrow(as.matrix(x))
  if (!(length(val) %in% c(1, n_mcmc))) {
    stop("val must have length 1 or number of mcmc samples.")
  }
  if (length(val) == 1) {
    val <- rep(val, n_mcmc)
  }
  val
}

#' @param x mcmc output
#' @param n_chains number of MCMC chains
#' @param times times in the data
assert_mcmc_format <- function(x, n_chains, times = NULL) {
  mcmc_names <- c("doses", "times", "w_prior", "w_post")
  attr_names <- c("names", "response_type", "doses", "model_index", "class")
  attr_names_mod <- c("class", "doses")
  if (!is.null(times)) {
    mcmc_names <- c("times", mcmc_names)
    attr_names <- c("longitudinal_model", attr_names)
    attr_names_mod <- c(attr_names_mod, "times", "longitudinal_model", "t_max")
  }
  expect_true(all(mcmc_names %in% names(x)))
  expect_equal(length(x$mod), n_chains)
  expect_true(all(attr_names %in% names(attributes(x))))
  expect_true(all(attr_names_mod %in% names(attributes(x$mod))))
}

#' @param mcmc output from dreamer_mcmc() with a single model named "mod"
#' @param dose a single dose at which to use in posterior()
#' @param b1,b2,b3,b4 vectors to replace MCMC samples
#' @param true_responses the true mean responses assuming b1 through b4
#'   (should be an expression, to be evaluated in test_posterior_impl,
#'   so it should be a function of objects in that environment)
test_posterior <- function(
  mcmc,
  prob,
  doses,
  times = NULL,
  reference_dose = NULL,
  b1 = NULL,
  b2 = NULL,
  b3 = NULL,
  b4 = NULL,
  a = NULL,
  c1 = NULL,
  c2 = NULL,
  d1 = NULL,
  d2 = NULL,
  gam = NULL,
  `b1[1]` = NULL, # nolint
  `b1[2]` = NULL, # nolint
  `b1[3]` = NULL, # nolint
  true_responses
) {
  mcmc <- mcmc %>%
    replace_mcmc("mod", "b1", b1) %>%
    replace_mcmc("mod", "b2", b2) %>%
    replace_mcmc("mod", "b3", b3) %>%
    replace_mcmc("mod", "b4", b4) %>%
    replace_mcmc("mod", "a", a) %>%
    replace_mcmc("mod", "c1", c1) %>%
    replace_mcmc("mod", "c2", c2) %>%
    replace_mcmc("mod", "d1", d1) %>%
    replace_mcmc("mod", "d2", d2) %>%
    replace_mcmc("mod", "gam", gam) %>%
    replace_mcmc("mod", "b1[1]", `b1[1]`) %>%
    replace_mcmc("mod", "b1[2]", `b1[2]`) %>%
    replace_mcmc("mod", "b1[3]", `b1[3]`)
  post <- posterior(
    mcmc,
    prob = prob,
    doses = doses,
    times = times,
    reference_dose = reference_dose,
    return_samples = TRUE
  )
  the_grid <- tidyr::expand_grid(dose = doses, time = times)
  purrr::pwalk(
    the_grid,
    test_posterior_impl,
    mcmc = mcmc,
    reference_dose = reference_dose,
    post = post,
    prob = prob,
    b1 = b1,
    b2 = b2,
    b3 = b3,
    b4 = b4,
    a = a,
    c1 = c1,
    c2 = c2,
    d1 = d1,
    d2 = d2,
    gam = gam,
    `b1[1]` = `b1[1]`,
    `b1[2]` = `b1[2]`,
    `b1[3]` = `b1[3]`,
    true_responses = true_responses
  )
}

test_posterior_impl <- function(
  mcmc,
  post,
  prob,
  dose,
  time = NULL,
  reference_dose = NULL,
  b1 = NULL,
  b2 = NULL,
  b3 = NULL,
  b4 = NULL,
  a = NULL,
  c1 = NULL,
  c2 = NULL,
  d1 = NULL,
  d2 = NULL,
  gam = NULL,
  `b1[1]` = NULL, # nolint 
  `b1[2]` = NULL, # nolint
  `b1[3]` = NULL, # nolint 
  true_responses
) {
  true_responses <- eval(true_responses)
  post <- subset_post(post, dose, time)
  true_mean <- mean(true_responses)
  true_q <- quantile(true_responses, prob, names = FALSE)
  lb <- sprintf("%.2f%%", 100 * prob[1])
  ub <- sprintf("%.2f%%", 100 * prob[2])
  expect_equal(post$samps$mean_response, true_responses)
  expect_equal(post$stats$mean, true_mean)
  expect_equal(c(post$stats[[lb]], post$stats[[ub]]), true_q)
}

subset_post <- function(post, dose, time) {
  post$stats <- dplyr::filter(post$stats, dose == !!dose)
  post$samps <- dplyr::filter(post$samps, dose == !!dose)
  if (!is.null(time)) {
    post$stats <- dplyr::filter(post$stats, time == !!time)
    post$samps <- dplyr::filter(post$samps, time == !!time)
  }
  post
}
