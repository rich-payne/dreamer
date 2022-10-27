assert_doses <- function(doses) {
  if (is.null(doses)) rlang::abort("Please specify doses.", class = "dreamer")
}

assert_data_reference_dose <- function(data, reference_dose) {
  if (!is.null(data) && !is.null(reference_dose)) {
    rlang::abort(
      "Plotting of placebo adjusted data not supported; omit 'data' argument.",
      class = "dreamer"
    )
  }
}

assert_no_dots <- function(function_name, ...) {
  args <- list(...)
  if (length(args) > 0) {
    msg <- paste0(
      "The following arguments are not supported for :\n",
      paste0(names(args), collapse = ", ")
    )
    rlang::abort(msg, class = "dreamer")
  }
}

assert_reference_dose <- function(reference_dose) {
  if (is.null(reference_dose)) return(NULL)
  if (length(reference_dose) != 1) {
    rlang::abort(
      "Argument `reference_dose` must have length 1 or be NULL.",
      class = "dreamer"
    )
  }
}

assert_model_names <- function(x) {
  for (i in seq_along(x)) {
    if (isTRUE(names(x)[i] == "") || is.null(names(x)[i])) {
      rlang::abort(
        "All dreamer model arguments must be named in dreamer_mcmc().",
        class = "dreamer"
      )
    }
  }
  if (any(duplicated(names(x)))) {
    rlang::abort("Duplicate model names are not allowed.", class = "dreamer")
  }
}

assert_names <- function(x, y, msg) {
  if (!all(names(x) == names(y))) {
    stop(msg, call. = FALSE)
  }
}

assert_w_prior <- function(w_prior) {
  w_total <- sum(w_prior)
  if (!isTRUE(all.equal(w_total, 1))) {
    stop(
      "Sum of w_prior for all models must add to 1: ", w_total,
      call. = FALSE
    )
  }
}
