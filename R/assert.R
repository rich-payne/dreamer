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
