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

#' @export
print.dreamer_continuous <- function(x, ...) {
  print_model(x, "continuous")
}

#' @export
print.dreamer_binary <- function(x, ...) {
  print_model(x, "binary")
}

print_model <- function(x, response_type) {
  hyparms <- x
  hyparms$longitudinal <- NULL
  hyparms <- get_hyparms(hyparms) %>%
    knitr::kable(format = "pipe")
  long <- get_hyparms(x$longitudinal)
  if (!is.null(long)) {
    long <- knitr::kable(long, format = "pipe")
  }
  cat("-------------\n")
  cat("dreamer Model\n")
  cat("-------------\n")
  cat(paste0("type    : ", attr(x, "type"), "\n"))
  cat(paste0("response: ", response_type, "\n"))
  cat("\n")
  cat("Dose Response\n")
  cat(hyparms, sep = "\n")
  
  if (!is.null(long)) {
    cat("\n")
    cat("Longitudinal\n")
    cat(long, sep = "\n")
  }
}

get_hyparms <- function(hyparms) {
  if (is.null(hyparms)) return()
  hyparms <- data.frame(
    hyperparameter = names(hyparms),
    value = unlist(hyparms, use.names = FALSE)
  )
}
