#' @title Traceplots
#' @description Produces traceplots for each parameter for each model.
#' @param x output from a call to `dreamer_mcmc()`.
#' @example man/examples/ex-plot_trace.R
#' @export
plot_trace <- function(x) {
  UseMethod("plot_trace", x)
}

#' @export
plot_trace.dreamer <- function(x) { #nolint
  graphics::par(mfrow = c(3, 2))
  coda::traceplot(x)
}

#' @export
plot_trace.dreamer_bma <- function(x) { #nolint
  ind <- vapply(
    x,
    function(model) any(grepl("mcmc.list", class(model))),
    logical(1)
  ) %>%
    which()
  graphics::par(mfrow = c(3, 2))
  for (i in ind) {
    for (j in seq_len(length(x[[i]]))) {
      colnames(x[[i]][[j]]) <- paste0(names(x)[i], "_", colnames(x[[i]][[j]]))
    }
    coda::traceplot(x[[i]])
  }
}
