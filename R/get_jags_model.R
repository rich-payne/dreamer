get_jags_model <- function(model, data) {
  mod_likelihood <- ""
  if (!is.null(data)) {
    mod_likelihood <- get_jags_like_file(model) %>% get_jags_code()
  }
  mod_prior <- get_jags_prior_file(model) %>% get_jags_code()
  long_mod <- 1
  long_prior <- ""
  if (!is.null(model$longitudinal)) {
    long_prior <- get_long_prior(model$longitudinal)
    long_mod <- get_long_model(model$longitudinal)
  }
  mod_likelihood <- gsub("LONGITUDINAL", long_mod, mod_likelihood) %>%
    set_link(model)
  create_jags_model(mod_likelihood, mod_prior, long_prior)
}

set_link <- function(mod_likelihood, model) {
  if (inherits(model, "dreamer_binary")) {
    mod_likelihood <- gsub("LINK", model$link, mod_likelihood)
  }
  return(mod_likelihood)
}

get_long_prior <- function(longitudinal) {
  long_prior <- get_longitudinal_prior_file(longitudinal) %>%
    get_jags_code()
  long_prior_intercept <- get_jags_code(
    "jags/prior_longitudinal_global_intercept.jags"
  )
  long_prior <- paste0(long_prior, long_prior_intercept, collapse = "\n")
  return(long_prior)
}

get_long_model <- function(longitudinal) {
  long_mod <- get_longitudinal_likelihood_file(longitudinal) %>%
    get_jags_code()
  return(long_mod)
}

create_jags_model <- function(...) {
  inputs <- list(...)
  inputs$sep <- "\n"
  model <- do.call(paste, inputs)
  paste("model{", model, "}", sep = "\n")
}

get_jags_code <- function(file) {
  system.file(
    file,
    package = "dreamer",
    mustWork = TRUE
  ) %>%
    readLines() %>%
    paste0(collapse = "\n")
}

get_jags_like_file <- function(x) {
  UseMethod("get_jags_like_file", x)
}

#' @export
get_jags_like_file.dreamer_linear <- function(x) { #nolint
  "jags/likelihood_linear.jags"
}

#' @export
get_jags_like_file.dreamer_quad <- function(x) { #nolint
  "jags/likelihood_quad.jags"
}

#' @export
get_jags_like_file.dreamer_loglinear <- function(x) { #nolint
  "jags/likelihood_loglinear.jags"
}

#' @export
get_jags_like_file.dreamer_logquad <- function(x) { #nolint
  "jags/likelihood_logquad.jags"
}

#' @export
get_jags_like_file.dreamer_emax <- function(x) { #nolint
  "jags/likelihood_emax.jags"
}

#' @export
get_jags_like_file.dreamer_exp <- function(x) {
  "jags/likelihood_exp.jags"
}

#' @export
get_jags_like_file.dreamer_beta <- function(x) { #nolint
  "jags/likelihood_beta.jags"
}

#' @export
get_jags_like_file.dreamer_independent_continuous <- function(x) { #nolint
  "jags/likelihood_independent.jags"
}

#' @export
get_jags_like_file.dreamer_linear_binary <- function(x) { #nolint
  "jags/likelihood_linear_binary.jags"
}

#' @export
get_jags_like_file.dreamer_quad_binary <- function(x) { #nolint
  "jags/likelihood_quad_binary.jags"
}

#' @export
get_jags_like_file.dreamer_loglinear_binary <- function(x) { #nolint
  "jags/likelihood_loglinear_binary.jags"
}

#' @export
get_jags_like_file.dreamer_logquad_binary <- function(x) { #nolint
  "jags/likelihood_logquad_binary.jags"
}

#' @export
get_jags_like_file.dreamer_emax_binary <- function(x) { #nolint
  "jags/likelihood_emax_binary.jags"
}

#' @export
get_jags_like_file.dreamer_exp_binary <- function(x) { #nolint
  "jags/likelihood_exp_binary.jags"
}

#' @export
get_jags_like_file.dreamer_beta_binary <- function(x) { #nolint
  "jags/likelihood_beta_binary.jags"
}

#' @export
get_jags_like_file.dreamer_independent_binary <- function(x) { #nolint
  "jags/likelihood_independent_binary.jags"
}

get_jags_prior_file <- function(x) {
  UseMethod("get_jags_prior_file", x)
}

#' @export
get_jags_prior_file.dreamer_linear <- function(x) { #nolint
  "jags/prior_linear.jags"
}

#' @export
get_jags_prior_file.dreamer_quad <- function(x) { #nolint
  "jags/prior_quad.jags"
}

#' @export
get_jags_prior_file.dreamer_loglinear <- function(x) { #nolint
  "jags/prior_loglinear.jags"
}

#' @export
get_jags_prior_file.dreamer_logquad <- function(x) { #nolint
  "jags/prior_logquad.jags"
}

#' @export
get_jags_prior_file.dreamer_emax <- function(x) { #nolint
  "jags/prior_emax.jags"
}

#' @export
get_jags_prior_file.dreamer_exp <- function(x) { #nolint
  "jags/prior_exp.jags"
}

#' @export
get_jags_prior_file.dreamer_beta <- function(x) { #nolint
  "jags/prior_beta.jags"
}

#' @export
get_jags_prior_file.dreamer_independent_continuous <- function(x) { #nolint
  "jags/prior_independent.jags"
}

#' @export
get_jags_prior_file.dreamer_linear_binary <- function(x) { #nolint
  "jags/prior_linear_binary.jags"
}

#' @export
get_jags_prior_file.dreamer_quad_binary <- function(x) { #nolint
  "jags/prior_quad_binary.jags"
}

#' @export
get_jags_prior_file.dreamer_loglinear_binary <- function(x) { #nolint
  "jags/prior_loglinear_binary.jags"
}

#' @export
get_jags_prior_file.dreamer_logquad_binary <- function(x) { #nolint
  "jags/prior_logquad_binary.jags"
}

#' @export
get_jags_prior_file.dreamer_emax_binary <- function(x) { #nolint
  "jags/prior_emax_binary.jags"
}

#' @export
get_jags_prior_file.dreamer_exp_binary <- function(x) { #nolint
  "jags/prior_exp_binary.jags"
}

#' @export
get_jags_prior_file.dreamer_beta_binary <- function(x) { #nolint
  "jags/prior_beta_binary.jags"
}

#' @export
get_jags_prior_file.dreamer_independent_binary <- function(x) { #nolint
  "jags/prior_independent_binary.jags"
}


get_vnames <- function(x) {
  c(get_vnames_model(x), get_longitudinal_vnames(x$longitudinal))
}

get_vnames_model <- function(x) {
  UseMethod("get_vnames_model", x)
}

#' @export
get_vnames_model.dreamer_linear <- function(x) { #nolint
  c("b1", "b2", "sigma")
}

#' @export
get_vnames_model.dreamer_quad <- function(x) {
  c("b1", "b2", "b3", "sigma")
}

#' @export
get_vnames_model.dreamer_loglinear <- function(x) { #nolint
  c("b1", "b2", "sigma")
}

#' @export
get_vnames_model.dreamer_logquad <- function(x) { #nolint
  c("b1", "b2", "b3", "sigma")
}

#' @export
get_vnames_model.dreamer_emax <- function(x) {
  c("b1", "b2", "b3", "b4", "sigma")
}

#' @export
get_vnames_model.dreamer_exp <- function(x) {
  c("b1", "b2", "b3", "sigma")
}

#' @export
get_vnames_model.dreamer_beta <- function(x) {
  c("b1", "b2", "b3", "b4", "sigma")
}

#' @export
get_vnames_model.dreamer_independent_continuous <- function(x) { #nolint
  c("b1", "sigma")
}

#' @export
get_vnames_model.dreamer_linear_binary <- function(x) { #nolint
  c("b1", "b2")
}

#' @export
get_vnames_model.dreamer_quad_binary <- function(x) { #nolint
  c("b1", "b2", "b3")
}

#' @export
get_vnames_model.dreamer_loglinear_binary <- function(x) { #nolint
  c("b1", "b2")
}

#' @export
get_vnames_model.dreamer_logquad_binary <- function(x) { #nolint
  c("b1", "b2", "b3")
}

#' @export
get_vnames_model.dreamer_emax_binary <- function(x) { #nolint
  c("b1", "b2", "b3", "b4")
}

#' @export
get_vnames_model.dreamer_exp_binary <- function(x) { #nolint
  c("b1", "b2", "b3")
}

#' @export
get_vnames_model.dreamer_beta_binary <- function(x) { #nolint
  c("b1", "b2", "b3", "b4")
}

#' @export
get_vnames_model.dreamer_independent_binary <- function(x) { #nolint
  c("b1")
}

get_longitudinal_likelihood_file <- function(x) {  #nolint
  UseMethod("get_longitudinal_likelihood_file", x)
}

#' @export
get_longitudinal_likelihood_file.dreamer_longitudinal_linear <- function(x) { #nolint
  "jags/likelihood_longitudinal_linear.jags"
}

#' @export
get_longitudinal_likelihood_file.dreamer_longitudinal_itp <- function(x) { #nolint
  "jags/likelihood_longitudinal_itp.jags"
}

#' @export
get_longitudinal_likelihood_file.dreamer_longitudinal_idp <- function(x) { #nolint
  "jags/likelihood_longitudinal_idp.jags"
}

get_longitudinal_prior_file <- function(x) { #nolint
  UseMethod("get_longitudinal_prior_file", x)
}

#' @export
get_longitudinal_prior_file.dreamer_longitudinal_linear <- function(x) { #nolint
  "jags/prior_longitudinal_linear.jags"
}

#' @export
get_longitudinal_prior_file.dreamer_longitudinal_itp <- function(x) { #nolint
  "jags/prior_longitudinal_itp.jags"
}

#' @export
get_longitudinal_prior_file.dreamer_longitudinal_idp <- function(x) { #nolint
  "jags/prior_longitudinal_idp.jags"
}

get_longitudinal_vnames <- function(x) {
  UseMethod("get_longitudinal_vnames", x)
}

#' @export
get_longitudinal_vnames.default <- function(x) NULL #nolint

#' @export
get_longitudinal_vnames.dreamer_longitudinal_linear <- function(x) { #nolint
  c("a")
}

#' @export
get_longitudinal_vnames.dreamer_longitudinal_itp <- function(x) { #nolint
  c("a", "c1")
}

#' @export
get_longitudinal_vnames.dreamer_longitudinal_idp <- function(x) { #nolint
  c("a", "c1", "c2", "d1", "d2", "gam")
}
