#' @title Calculate MCMC Diagnostics for Parameters
#' @description Calculate MCMC diagnostics for individual parameters.
#' @param x MCMC output from a dreamer model.
#' @example man/examples/ex-diagnostics.R
#' @export
diagnostics <- function(x) {
  UseMethod("diagnostics", x)
}

#' @export
diagnostics.dreamer_bma <- function(x) { #nolint
  model_index <- vapply(x, function(x) "mcmc.list" %in% class(x), logical(1))
  model_list <- x[model_index]
  diags <- purrr::map(
    model_list,
    function(y) {
      diagnostics(y)
    }
  )
  diags <- purrr::map2(
    diags,
    names(diags),
    function(x, y) dplyr::mutate(x, model = y)
  ) %>%
    do_call(getExportedValue("dplyr", "bind_rows")) %>%
    dplyr::select(.data$model, everything())
  return(diags)
}

#' @export
diagnostics.dreamer <- function(x) {
  diags <- coda::gelman.diag(x, multivariate = FALSE)$psrf %>%
    dplyr::as_tibble(rownames = "param") %>%
    dplyr::select(
      .data$param,
      gelman_point = .data$`Point est.`,
      gelman_upper = .data$`Upper C.I.`
    ) %>%
    dplyr::mutate(effective_size = coda::effectiveSize(x)[.data$param])
  return(diags)
}
