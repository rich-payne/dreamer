#' @import coda
#' @importFrom dplyr
#'   .data
#'   %>%
#'   across
#'   any_of
#'   desc
#'   everything
#'   n
#' @importFrom ggplot2
#'   aes
#'   aes_string
#'   element_text
#'   geom_errorbar
#'   geom_line
#'   geom_point
#'   ggplot
#'   labs
#'   scale_x_continuous
#'   theme
#' @importFrom rlang :=
#' @importFrom stats
#'   pnorm
#'   quantile
#'   rbinom
#'   rnorm
#'   dbinom
#'   dnorm
#'   var
#'   qnorm
utils::globalVariables(c(".")) # remove binding notes
