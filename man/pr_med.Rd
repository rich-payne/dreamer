% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pr_med.R
\name{pr_med}
\alias{pr_med}
\title{Pr(minimum efficacious dose)}
\usage{
pr_med(
  x,
  doses = attr(x, "doses"),
  csd = NULL,
  reference_dose = NULL,
  greater = TRUE,
  time = NULL
)
}
\arguments{
\item{x}{output from a call to \code{dreamer_mcmc()}.}

\item{doses}{the doses for which pr(MED) is to be calculated.}

\item{csd}{the treatment effect that is clinically relevant.}

\item{reference_dose}{a single dose that is used as the reference when
defining the MED relative to a dose (rather than in absolute terms).  When
\code{reference_dose} is specified, this function calculates the
posterior probability
that each dose is the smallest dose such that
(effect_dose - effect_reference_dose > CSD).}

\item{greater}{if \code{TRUE}, higher responses indicate better efficacy.  If
\code{FALSE}, lower responses indicate better efficacy.`}

\item{time}{the time (scalar) at which the Pr(MED) should be calculated.
Applies only to longitudinal models.}
}
\value{
A tibble listing each dose and the posterior probability that
each dose is the minimum efficacious dose.
}
\description{
Calculates the posterior probability that each specified doses
are the minimum effective dose in the set; i.e. the smallest
dose that has a clinically significant difference (CSD).
}
\examples{
set.seed(888)
data <- dreamer_data_linear(
  n_cohorts = c(20, 20, 20),
  dose = c(0, 3, 10),
  b1 = 1,
  b2 = 3,
  sigma = 5
)

# Bayesian model averaging
output <- dreamer_mcmc(
 data = data,
 n_adapt = 1e3,
 n_burn = 1e3,
 n_iter = 1e3,
 n_chains = 2,
 silent = FALSE,
 mod_linear = model_linear(
   mu_b1 = 0,
   sigma_b1 = 1,
   mu_b2 = 0,
   sigma_b2 = 1,
   shape = 1,
   rate = .001,
   w_prior = 1 / 2
 ),
 mod_quad = model_quad(
   mu_b1 = 0,
   sigma_b1 = 1,
   mu_b2 = 0,
   sigma_b2 = 1,
   mu_b3 = 0,
   sigma_b3 = 1,
   shape = 1,
   rate = .001,
   w_prior = 1 / 2
 )
)

pr_med(output, csd = 10)

# difference of two doses
pr_med(output, csd = 3, reference_dose = 0)

# single model
pr_med(output$mod_quad, csd = 10)
}
