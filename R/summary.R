#' @title Summarize Bayesian Model Averaging MCMC Output
#' @description Summarize parameter inference and convergence diagnostics.
#' @param object a dreamer MCMC object.
#' @param ... additional arguments (which are ignored).
#' @return Returns a named list with elements `model_weights` and `summary`
#'   containing the prior and posterior weights for each model and inference
#'   on parameters for each model as well as MCMC diagnostics.
#' @example man/examples/ex-summary.R
#' @export
summary.dreamer_bma <- function(object, ...) {
  assert_no_dots("summary.dreamer_bma", ...)
  model_index <- vapply(
    object,
    function(x) "mcmc.list" %in% class(x),
    logical(1)
  )
  model_list <- object[model_index]
  w_post_df <- dplyr::tibble(
    model = names(object$w_post),
    w_post = object$w_post,
    posterior_weight = sprintf("%.1f%%", object$w_post * 100),
    prior_weight = sprintf("%.1f%%", object$w_prior * 100)
  ) %>%
    dplyr::arrange(desc(.data$w_post)) %>%
    dplyr::select(- .data$w_post)
  sumry <- lapply(model_list, summary)
  sumry <- purrr::map2(
    sumry,
    names(sumry),
    function(x, y) dplyr::mutate(x, model = y)
  ) %>%
    do_call(getExportedValue("dplyr", "bind_rows")) %>%
    dplyr::select(.data$model, everything())
  the_summary <- list(
    model_weights = w_post_df,
    summary = sumry
  )
  return(the_summary)
}

#' @title Summarize Model Output
#' @description Produces summaries for inference and diagnosing MCMC chains.
#' @param object MCMC output from a dreamer model.
#' @param ... additional arguments which are ignored.
#' @return A tibble with inference and diagnostics information for each
#'   parameter.
#' @example man/examples/ex-summary.R
#' @export
summary.dreamer_mcmc <- function(object, ...) {
  assert_no_dots("summary.dreamer_mcmc", ...)
  diags <- diagnostics(object)
  out <- NextMethod("summary", object) # use coda summary
  param <- rownames(out$statistics)
  if (is.null(param)) { # if only one parameter to be estimated
    param <- colnames(object[[1]])[1]
  }
  if (is.null(dim(out$statistics))) {
    cnames <- names(out$statistics)
    dim(out$statistics) <- c(1, length(out$statistics))
    colnames(out$statistics) <- cnames
  }
  if (is.null(dim(out$quantiles))) {
    cnames <- names(out$quantiles)
    dim(out$quantiles) <- c(1, length(out$quantiles))
    colnames(out$quantiles) <- cnames
  }
  inference <- cbind(out$statistics, out$quantiles) %>%
    dplyr::as_tibble() %>%
    dplyr::bind_cols(param = param) %>%
    dplyr::select(
      param,
      mean = .data$Mean,
      sd = .data$SD,
      se = .data$`Naive SE`,
      se_ts = .data$`Time-series SE`,
      dplyr::everything()
    )
  output <- dplyr::full_join(inference, diags, by = "param")
  return(output)
}

do_call <- function(args, what) {
  do.call(what = what, args = args)
}
