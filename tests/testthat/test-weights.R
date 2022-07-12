t_max <- 5
long_mods <- tibble::tibble(
  longitudinal = list(
    NULL,
    rlang::expr(
      model_longitudinal_linear(
        mu_a = 0,
        sigma_a = 1,
        t_max = !!t_max
      )
    ),
    rlang::expr(
      model_longitudinal_itp(
        mu_a = 0,
        sigma_a = 1,
        a_c1 = 0,
        b_c1 = 1,
        t_max = !!t_max
      )
    ),
    rlang::expr(
      model_longitudinal_idp(
        mu_a = 0,
        sigma_a = 1,
        a_c1 = 0,
        b_c1 = 1,
        a_c2 = - 1,
        b_c2 = 0,
        t_max = !!t_max
      )
    )
  )
)

parms <- list(
  a = .5,
  b1 = .3,
  b2 = .4,
  b3 = .5,
  b4 = .6,
  sigma = 4,
  c1 = .25,
  c2 = - .5,
  d1 = 1.1,
  d2 = 3.5,
  gam = .75
)

replace_params <- function(x, model_name, ...) {
  vars <- list(...)
  varnames <- names(vars)
  mcmc <- x[[model_name]][[1]]
  for (i in seq_len(length(vars))) {
    if (varnames[i] %in% colnames(mcmc)) {
      x[[model_name]][[1]][, varnames[i]] <- vars[[i]]
      mcmc <- x[[model_name]][[1]]
    } else {
      newmcmc <- cbind(
        mcmc,
        matrix(
          vars[[i]],
          nrow = nrow(mcmc),
          dimnames = list(NULL, varnames[i])
        )
      )
      attr(newmcmc, "mcpar") <- attr(mcmc, "mcpar")
      class(newmcmc) <- class(mcmc)
      x[[model_name]][[1]] <- newmcmc
      mcmc <- newmcmc
    }
  }
  return(x)
}

test_that("model weights", {
  set.seed(8888111)

  w_prior <- 1 / 7
  dose_mods <- tibble::tibble(
    mod_lin = list(rlang::expr(
      model_linear(
        mu_b1 = 0,
        sigma_b1 = 1,
        mu_b2 = 0,
        sigma_b2 = 1,
        shape = 1,
        rate = .01,
        w_prior = !!w_prior
      )
    )),
    mod_quad = list(rlang::expr(
      model_quad(
        mu_b1 = 0,
        sigma_b1 = 1,
        mu_b2 = 0,
        sigma_b2 = 1,
        mu_b3 = 0,
        sigma_b3 = 1,
        shape = 1,
        rate = .01,
        w_prior = !!w_prior
      )
    )),
    mod_loglin = list(rlang::expr(
      model_loglinear(
        mu_b1 = 0,
        sigma_b1 = 1,
        mu_b2 = 0,
        sigma_b2 = 1,
        shape = 1,
        rate = .01,
        w_prior = !!w_prior
      )
    )),
    mod_logquad = list(rlang::expr(
      model_logquad(
        mu_b1 = 0,
        sigma_b1 = 1,
        mu_b2 = 0,
        sigma_b2 = 1,
        mu_b3 = 0,
        sigma_b3 = 1,
        shape = 1,
        rate = .01,
        w_prior = !!w_prior
      )
    )),
    mod_emax = list(rlang::expr(
      model_emax(
        mu_b1 = 0,
        sigma_b1 = 1,
        mu_b2 = 0,
        sigma_b2 = 1,
        mu_b3 = 0,
        sigma_b3 = 1,
        mu_b4 = 0,
        sigma_b4 = 1,
        shape = 1,
        rate = .01,
        w_prior = !!w_prior
      )
    )),
    mod_exp = list(rlang::expr(
      model_exp(
        mu_b1 = 0,
        sigma_b1 = 1,
        mu_b2 = 0,
        sigma_b2 = 1,
        mu_b3 = 0,
        sigma_b3 = 1,
        shape = 1,
        rate = .01,
        w_prior = !!w_prior
      )
    )),
    mod_beta = list(rlang::expr(
      model_beta(
        mu_b1 = 0,
        sigma_b1 = 1,
        mu_b2 = 0,
        sigma_b2 = 1,
        mu_b3 = 0,
        sigma_b3 = 1,
        mu_b4 = 0,
        sigma_b4 = 1,
        shape = 1,
        rate = .01,
        w_prior = !!w_prior
      )
    ))
  ) %>%
    tidyr::pivot_longer(everything(), names_to = "model")

  scens <- tidyr::expand_grid(dose_mods, long_mods) %>%
    dplyr::mutate(
      value = purrr::map2(
        value,
        longitudinal,
        ~ eval(rlang::call_modify(.x, longitudinal = .y))
      )
    ) %>%
    tidyr::pivot_wider(names_from = "model", values_from = "value")

  check_weights <- function(data, parms, ...) {
    out <- rlang::exec(
      dreamer_mcmc,
      data = data,
      n_iter = 2,
      n_chains = 1,
      silent = TRUE,
      convergence_warn = FALSE,
      ...
    )

    n <- nrow(data)
    response <- data$response
    dose <- data$dose
    a <- parms$a
    b1 <- parms$b1
    b2 <- parms$b2
    b3 <- parms$b3
    b4 <- parms$b4
    sigma <- parms$sigma
    c1 <- parms$c1
    c2 <- parms$c2
    d1 <- parms$d1
    d2 <- parms$d2
    gam <- parms$gam
    scale <- attr(out$mod_beta, "scale")
    t_max <- attr(out$mod_lin, "t_max")

    lmod <- attr(out$mod_lin, "longitudinal")
    if (is.null(lmod)) {
      time_perc <- 1
      a <- 0
      nlp <- 0
    } else if (any(lmod == "dreamer_longitudinal_linear")) {
      time_perc <- data$time / t_max
      nlp <- 1
    } else if (any(lmod == "dreamer_longitudinal_itp")) {
      time_perc <- ((1 - exp(- c1 * data$time)) / (1 - exp(- c1 * t_max)))
      nlp <- 2
    } else if (any(lmod == "dreamer_longitudinal_idp")) {
      time_perc <- (((1 - exp(- c1 * data$time)) / (1 - exp(- c1 * d1))) *
        (data$time < d1) + (1 - gam * ((1 - exp(- c2 * (data$time - d1))) /
        (1 - exp(- c2 * (d2 - d1))))) *
        ((data$time >= d1) & (data$time <= d2)) + (1 - gam) * (data$time > d2))
      nlp <- 6
    }

    wl_lin <- sum(
      dnorm(
        response,
        a + (b1 + b2 * dose) * time_perc,
        sigma,
        log = TRUE
      )
    ) - (3 + nlp) / 2
    wl_quad <- sum(
      dnorm(
        response,
        a + (b1 + b2 * dose + b3 * dose^2) * time_perc,
        sigma,
        log = TRUE
      )
    ) - (4 + nlp) / 2
    wl_loglin <- sum(
      dnorm(
        response,
        a + (b1 + b2 * log(dose + 1)) * time_perc,
        sigma,
        log = TRUE
      )
    ) - (3 + nlp) / 2
    wl_logquad <- sum(
      dnorm(
        response,
        a + (b1 + b2 * log(dose + 1) + b3 * log(dose + 1)^2) * time_perc,
        sigma,
        log = TRUE
      )
    ) - (4 + nlp) / 2
    wl_emax <- sum(
      dnorm(
        response,
        a + (b1 + (b2 - b1) * dose ^ b4 / (exp(b3 * b4) + dose ^ b4)) *
          time_perc,
        sigma,
        log = TRUE
      )
    ) - (5 + nlp) / 2
    wl_exp <- sum(
      dnorm(
        response,
        a + (b1 + b2 * (1 - exp(- b3 * dose))) * time_perc,
        sigma,
        log = TRUE
      )
    ) - (4 + nlp) / 2
    wl_beta <- sum(
      dnorm(
        response,
        a + (b1 + b2 * ((b3 + b4) ^ (b3 + b4)) / (b3 ^ b3 * b4 ^ b4) *
          (dose / scale) ^ b3 * (1 - dose / scale) ^ b4) * time_perc,
        sigma,
        log = TRUE
      )
    ) - (5 + nlp) / 2

    wl <- log(w_prior) +
      c(wl_lin, wl_quad, wl_loglin, wl_logquad, wl_emax, wl_exp, wl_beta)
    the_max <- max(wl)
    denom <- the_max + log(sum(exp(wl - the_max)))
    true_weights <- exp(wl - denom)

    out2 <- out %>%
      {rlang::exec("replace_params", ., model_name = "mod_lin", !!!parms)} %>% #nolint
      {rlang::exec("replace_params", ., model_name = "mod_quad", !!!parms)} %>% #nolint
      {rlang::exec("replace_params", ., model_name = "mod_loglin", !!!parms)} %>% #nolint
      {rlang::exec("replace_params", ., model_name = "mod_logquad", !!!parms)} %>% #nolint
      {rlang::exec("replace_params", ., model_name = "mod_emax", !!!parms)} %>% #nolint
      {rlang::exec("replace_params", ., model_name = "mod_exp", !!!parms)} %>% #nolint
      {rlang::exec("replace_params", ., model_name = "mod_beta", !!!parms)}  #nolint

    mod_index <- vapply(
      out2,
      function(xx) "dreamer_mcmc" %in% class(xx),
      logical(1)
    )
    w_priors <- rep(w_prior, 7)
    names(w_priors) <- names(out2)[mod_index]
    names(true_weights) <- names(out2)[mod_index]
    dreamer_weights <- dreamer:::dreamer_post_weights(
      out2[mod_index],
      w_priors,
      data = data
    )
    expect_equal(dreamer_weights$w_post, true_weights)
  }

  set.seed(885)
  data <- dreamer_data_linear(
    n_cohorts = c(10, 15, 20, 25, 30),
    dose = c(1:5),
    b1 = .5,
    b2 = .3,
    sigma = .5,
    times = 1:5,
    longitudinal = "linear",
    a = .5,
    t_max = t_max
  )

  out <- purrr::pmap(
    dplyr::select(scens, - longitudinal),
    check_weights,
    parms = parms,
    data = data
  )
})

test_that("model weights (binary)", {
  set.seed(8888)
  w_prior <- 1 / 7
  dose_mods_bin <- tibble::tibble(
    mod_lin = list(rlang::expr(
      model_linear_binary(
        mu_b1 = 0,
        sigma_b1 = 1,
        mu_b2 = 0,
        sigma_b2 = 1,
        w_prior = !!w_prior
      )
    )),
    mod_quad = list(rlang::expr(
      model_quad_binary(
        mu_b1 = 0,
        sigma_b1 = 1,
        mu_b2 = 0,
        sigma_b2 = 1,
        mu_b3 = 0,
        sigma_b3 = 1,
        w_prior = !!w_prior
      )
    )),
    mod_loglin = list(rlang::expr(
      model_loglinear_binary(
        mu_b1 = 0,
        sigma_b1 = 1,
        mu_b2 = 0,
        sigma_b2 = 1,
        w_prior = !!w_prior
      )
    )),
    mod_logquad = list(rlang::expr(
      model_logquad_binary(
        mu_b1 = 0,
        sigma_b1 = 1,
        mu_b2 = 0,
        sigma_b2 = 1,
        mu_b3 = 0,
        sigma_b3 = 1,
        w_prior = !!w_prior
      )
    )),
    mod_emax = list(rlang::expr(
      model_emax_binary(
        mu_b1 = 0,
        sigma_b1 = 1,
        mu_b2 = 0,
        sigma_b2 = 1,
        mu_b3 = 0,
        sigma_b3 = 1,
        mu_b4 = 0,
        sigma_b4 = 1,
        w_prior = !!w_prior
      )
    )),
    mod_exp = list(rlang::expr(
      model_exp_binary(
        mu_b1 = 0,
        sigma_b1 = 1,
        mu_b2 = 0,
        sigma_b2 = 1,
        mu_b3 = 0,
        sigma_b3 = 1,
        w_prior = !!w_prior
      )
    )),
    mod_beta = list(rlang::expr(
      model_beta_binary(
        mu_b1 = 0,
        sigma_b1 = 1,
        mu_b2 = 0,
        sigma_b2 = 1,
        mu_b3 = 0,
        sigma_b3 = 1,
        mu_b4 = 0,
        sigma_b4 = 1,
        w_prior = !!w_prior
      )
    ))
  ) %>%
    tidyr::pivot_longer(everything(), names_to = "model") %>%
    tidyr::expand_grid(link = c("logit", "probit"))

  scens_bin <- tidyr::expand_grid(dose_mods_bin, long_mods) %>%
    dplyr::mutate(
      value = purrr::pmap(
        .,
        function(value, link, longitudinal, ...) {
          eval(
            rlang::call_modify(value, longitudinal = longitudinal, link = link)
          )
        }
      )
    ) %>%
    tidyr::pivot_wider(names_from = "model", values_from = "value")

  set.seed(885)
  data_bin <- dreamer_data_linear_binary(
    n_cohorts = c(10, 15, 20, 25, 30),
    dose = c(1:5),
    b1 = - 3,
    b2 = .25,
    times = 1:5,
    longitudinal = "linear",
    a = - .1,
    t_max = t_max,
    link = "logit"
  ) %>%
    dplyr::mutate(n = 1)

  check_weights_binary <- function(data, parms, ...) {
    out <- rlang::exec(
      dreamer_mcmc,
      data = data,
      n_iter = 2,
      n_chains = 1,
      silent = TRUE,
      convergence_warn = FALSE,
      ...
    )

    n <- nrow(data)
    response <- data$response
    dose <- data$dose
    a <- parms$a
    b1 <- parms$b1
    b2 <- parms$b2
    b3 <- parms$b3
    b4 <- parms$b4
    sigma <- parms$sigma
    c1 <- parms$c1
    c2 <- parms$c2
    d1 <- parms$d1
    d2 <- parms$d2
    gam <- parms$gam
    scale <- attr(out$mod_beta, "scale")
    t_max <- attr(out$mod_lin, "t_max")
    link <- attr(out$mod_lin, "link")
    ilink <- getFromNamespace(paste0("i", link), "dreamer")

    lmod <- attr(out$mod_lin, "longitudinal")
    if (is.null(lmod)) {
      time_perc <- 1
      a <- 0
      nlp <- 0
    } else if (any(lmod == "dreamer_longitudinal_linear")) {
      time_perc <- data$time / t_max
      nlp <- 1
    } else if (any(lmod == "dreamer_longitudinal_itp")) {
      time_perc <- ((1 - exp(- c1 * data$time)) / (1 - exp(- c1 * t_max)))
      nlp <- 2
    } else if (any(lmod == "dreamer_longitudinal_idp")) {
      time_perc <- (((1 - exp(- c1 * data$time)) / (1 - exp(- c1 * d1))) *
        (data$time < d1) + (1 - gam * ((1 - exp(- c2 * (data$time - d1))) /
        (1 - exp(- c2 * (d2 - d1))))) *
        ((data$time >= d1) & (data$time <= d2)) + (1 - gam) * (data$time > d2))
      nlp <- 6
    }

    wl_lin <- sum(
      dbinom(
        response,
        1,
        ilink(a + (b1 + b2 * dose) * time_perc),
        log = TRUE
      )
    ) - (2 + nlp) / 2
    wl_quad <- sum(
      dbinom(
        response,
        1,
        ilink(a + (b1 + b2 * dose + b3 * dose^2) * time_perc),
        log = TRUE
      )
    ) - (3 + nlp) / 2
    wl_loglin <- sum(
      dbinom(
        response,
        1,
        ilink(a + (b1 + b2 * log(dose + 1)) * time_perc),
        log = TRUE
      )
    ) - (2 + nlp) / 2
    wl_logquad <- sum(
      dbinom(
        response,
        1,
        ilink(a + (b1 + b2 * log(dose + 1) + b3 * log(dose + 1)^2) * time_perc),
        log = TRUE
      )
    ) - (3 + nlp) / 2
    wl_emax <- sum(
      dbinom(
        response,
        1,
        ilink(
          a + (b1 + (b2 - b1) * dose ^ b4 / (exp(b3 * b4) + dose ^ b4)) *
            time_perc
        ),
        log = TRUE
      )
    ) - (4 + nlp) / 2
    wl_exp <- sum(
      dbinom(
        response,
        1,
        ilink(a + (b1 + b2 * (1 - exp(- b3 * dose))) * time_perc),
        log = TRUE
      )
    ) - (3 + nlp) / 2
    wl_beta <- sum(
      dbinom(
        response,
        1,
        ilink(
          a + (b1 + b2 * ((b3 + b4) ^ (b3 + b4)) /
            (b3 ^ b3 * b4 ^ b4) *
            (dose / scale) ^ b3 * (1 - dose / scale) ^ b4) * time_perc
        ),
        log = TRUE)
      ) - (4 + nlp) / 2

    wl <- log(w_prior) +
      c(wl_lin, wl_quad, wl_loglin, wl_logquad, wl_emax, wl_exp, wl_beta)
    the_max <- max(wl)
    denom <- the_max + log(sum(exp(wl - the_max)))
    true_weights <- exp(wl - denom)

    out2 <- out %>%
      {rlang::exec("replace_params", ., model_name = "mod_lin", !!!parms)} %>% #nolint
      {rlang::exec("replace_params", ., model_name = "mod_quad", !!!parms)} %>% #nolint
      {rlang::exec("replace_params", ., model_name = "mod_loglin", !!!parms)} %>% #nolint
      {rlang::exec("replace_params", ., model_name = "mod_logquad", !!!parms)} %>% #nolint
      {rlang::exec("replace_params", ., model_name = "mod_emax", !!!parms)} %>% #nolint
      {rlang::exec("replace_params", ., model_name = "mod_exp", !!!parms)} %>% #nolint
      {rlang::exec("replace_params", ., model_name = "mod_beta", !!!parms)} #nolint

    mod_index <- vapply(
      out2,
      function(xx) "dreamer_mcmc_binary" %in% class(xx),
      logical(1)
    )
    w_priors <- rep(w_prior, 7)
    names(w_priors) <- names(out2)[mod_index]
    names(true_weights) <- names(out2)[mod_index]
    dreamer_weights <- dreamer:::dreamer_post_weights(
      out2[mod_index],
      w_priors,
      data = data
    )
    expect_equal(dreamer_weights$w_post, true_weights)
  }

  out_bin <- purrr::pmap(
    dplyr::select(scens_bin, - longitudinal, - link),
    check_weights_binary,
    parms = parms,
    data = data_bin
  )
})
