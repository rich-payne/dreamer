#' @export
print.dreamer_bma <- function(x, ...) {
  ellipsis::check_dots_empty()
  NextMethod("print", x)
  # cat("-----------------\n")
  # cat("dreamer Model Fit\n")
  # cat("-----------------\n")
  # cat(paste0("doses: ", paste(attr(x, "doses"), collapse = ", "), "\n"))
  # if (!is.null(attr(x, "times"))) {
  #   cat(paste0("times: ", paste0(attr(x, "times"), collapse = ", "), "\n"))
  # }
  weights <- dplyr::bind_cols(
    model = names(x$w_prior),
    `prior weight` = x$w_prior,
    `posterior weight` = x$w_post
  ) %>%
    knitr::kable(digits = 3, format = "pipe")
  cat("\n")
  cat(weights, sep = "\n")
  return(invisible())
}

#' @export
print.dreamer_mcmc <- function(x, ...) {
  ellipsis::check_dots_empty()
  cat("-----------------\n")
  cat("dreamer Model Fit\n")
  cat("-----------------\n")
  cat(paste0("doses: ", paste(attr(x, "doses"), collapse = ", "), "\n"))
  if (!is.null(attr(x, "times"))) {
    cat(paste0("times: ", paste0(attr(x, "times"), collapse = ", "), "\n"))
  }
  post <- knitr::kable(posterior(x)$stats, digits = 3, format = "pipe")
  cat("\n")
  cat(post, sep = "\n")
  return(invisible())
}