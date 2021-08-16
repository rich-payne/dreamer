get_jags_data <- function(model, data, doses) {
  prior_hyperparms <- get_jags_data_hyperparms(model, data) %>%
    add_long_hyperparms(model, data)
  jags_data <- prep_jags_data(model, data, doses = doses) %>%
    set_a_jags(model, data)
  out <- c(prior_hyperparms, jags_data)
  return(out)
}

add_long_hyperparms <- function(x, model, data) {
  out <- c(x, model$longitudinal)
  if (is.null(data)) {
    out <- remove_t_max(out, model$longitudinal)
  }
  return(out)
}

set_a_jags <- function(jags_data, model, data) {
  if (is.null(data)) {
    return(jags_data)
  }
  if (is.null(model$longitudinal)) {
    jags_data$a <- 0
  }
  return(jags_data)
}

get_jags_data_hyperparms <- function(model, ...) {
  UseMethod("get_jags_data_hyperparms", model)
}

get_jags_data_hyperparms.default <- function(model, data) { #nolint
  model[
    c("w_prior", "jags_rng", "jags_seed", "longitudinal", "link", "scale")
  ] <- NULL
  if (is.null(data)) {
    model <- model %>%
      remove_scale()
  }
  unclass(model)
}

get_jags_data_hyperparms.dreamer_independent <- function(model, data) { #nolint
  doses <- model$doses
  check_independent_model(data, doses)
  u_doses_data <- unique(data$dose)
  n_doses_data <- length(u_doses_data)
  if (is.null(doses)) {
    model$mu_b1 <- rep(model$mu_b1, n_doses_data)
    model$sigma_b1 <- rep(model$sigma_b1, n_doses_data)
  } else {
    d_ind <- order(unique(doses))
    doses <- doses[d_ind]
    check_doses_data(doses, data, u_doses_data)
    model$mu_b1 <- model$mu_b1[d_ind]
    model$sigma_b1 <- model$sigma_b1[d_ind]
  }
  model[c(
    "w_prior",
    "jags_rng",
    "jags_seed",
    "longitudinal",
    "doses",
    "link",
    "scale"
  )] <- NULL
  unclass(model)
}

check_doses_data <- function(doses, data, u_doses_data) {
  if (!is.null(data) && !isTRUE(all.equal(doses, sort(u_doses_data)))) {
    stop("Doses specified do not match doses in data", call. = FALSE)
  }
}

set_scale <- function(jags_data, model, data) {
  if (is.null(data)) return(jags_data)
  jags_data$scale <- get_scale(model, data)
  return(jags_data)
}

prep_jags_data <- function(model, data, ...) {
  UseMethod("prep_jags_data", model)
}

prep_jags_data.dreamer_binary <- function(model, data, ...) {
  prep_binary_jags_data(model, data)
}

prep_jags_data.dreamer_beta <- function(model, data, ...) {
  out <- prep_cont_jags_data(model, data)
  out <- out %>% set_scale(model, data)
  return(out)
}

prep_jags_data.dreamer_beta_binary <- function(model, data, ...) { #nolint
  out <- prep_binary_jags_data(model, data)
  out <- out %>% set_scale(model, data)
  return(out)
}

prep_binary_jags_data <- function(model, data) {
  if (is.null(data)) return(data)
  is_longitudinal <- !is.null(model$longitudinal)
  data <- aggregate_binary(data, is_longitudinal)
  jags_data <- list(
    y = data$response,
    dose = data$dose,
    n = data$n,
    n_obs = nrow(data)
  )
  if (is_longitudinal) {
    jags_data$time <- dplyr::pull(data, .data$time)
  }
  return(jags_data)
}

prep_jags_data.dreamer_continuous <- function(model, data, ...) { #nolint
  prep_cont_jags_data(model, data)
}

prep_jags_data.dreamer_independent_continuous <- function(model, data, doses) { #nolint
  if (is.null(data)) {
    jags_data <- list(n_doses = length(model$doses))
    return(jags_data)
  }
  jags_data <- prep_cont_jags_data(model, data)
  jags_data$dose_sufficient_index <- vapply(
    jags_data$dose_sufficient,
    function(yy) which(yy == doses),
    integer(1)
  )
  jags_data$n_doses <- length(unique(data$dose))
  jags_data$dose_sufficient <- NULL
  jags_data$doses <- NULL
  return(jags_data)
}

prep_jags_data.dreamer_independent_binary <- function(model, data, doses) { #nolint
  if (is.null(data)) {
    jags_data <- list(n_doses = length(model$doses))
    return(jags_data)
  }
  jags_data <- prep_binary_jags_data(model, data)
  jags_data$dose_index <- vapply(
    jags_data$dose,
    function(yy) which(yy == doses),
    integer(1)
  )
  jags_data$n_doses <- length(unique(data$dose))
  jags_data$dose <- NULL
  return(jags_data)
}

prep_cont_jags_data <- function(model, data) {
  if (is.null(data)) return(data)
  check_colnames(data)
  if (all(c("n", "sample_var") %in% colnames(data))) {
    time_var <- get_time_var(model)
    data_sufficient <- dplyr::select(
      data,
      .data$dose,
      ybar = .data$response,
      .data$sample_var,
      .data$n,
      any_of(!!time_var)
    )
  } else {
    data$n <- 1
    if (!is.null(model$longitudinal)) {
      data <- dplyr::group_by(data, .data$time)
    }
    data_sufficient <- dplyr::group_by(data, .data$dose, .add = TRUE) %>%
      dplyr::summarise(
        ybar = mean(.data$response),
        sample_var = var(.data$response),
        n = n()
      )
  }
  data_sufficient_vars <- dplyr::filter(data_sufficient, n > 1)

  jags_data <- list(
    n_means = data_sufficient$n,
    n_vars = data_sufficient_vars$n,
    ybar = data_sufficient$ybar,
    sample_var = data_sufficient_vars$sample_var,
    n_sufficient_means = nrow(data_sufficient),
    n_sufficient_vars = nrow(data_sufficient_vars),
    dose_sufficient = data_sufficient$dose
  )
  if (tibble::has_name(data_sufficient, "time")) {
    jags_data$time <- dplyr::pull(data_sufficient, .data$time)
  }
  return(jags_data)
}

get_time_var <- function(model) {
  if (is.null(model$longitudinal)) {
    return(NULL)
  } else {
    return("time")
  }
}

check_colnames <- function(data) {
  included <- c("n", "sample_var") %in% colnames(data)
  if (any(included) && !all(included)) {
    stop(
      "Please supply both 'n' and 'sample_var', or neither (raw data)",
      call. = FALSE
    )
  }
}

remove_t_max <- function(x, long_model) {
  UseMethod("remove_t_max", long_model)
}

#' @export
#' @keywords internal
remove_t_max.default <- function(x, ...) {
  x["t_max"] <- NULL
  return(x)
}

#' @export
#' @keywords internal
remove_t_max.dreamer_longitudinal_idp <- function(x, ...) { #nolint
  return(x)
}

remove_scale <- function(x) {
  x["scale"] <- NULL
  return(x)
}

aggregate_binary <- function(data, is_longitudinal) {
  if (is.null(data)) return(data)
  check_binary(data)
  if (!rlang::has_name(data, "n")) {
    data <- dplyr::mutate(data, n = 1)
  }
  if (is_longitudinal) {
    data <- dplyr::group_by(data, .data$time)
  }
  data <- dplyr::group_by(data, .data$dose, .add = TRUE) %>%
    dplyr::summarise(
      response = sum(.data$response),
      n = sum(.data$n)
    )
  return(data)
}

check_independent_model <- function(data, doses) {
  if (is.null(data) & is.null(doses)) {
    stop("specify data, or supply doses to independent model", call. = FALSE)
  }
}
