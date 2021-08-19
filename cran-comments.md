# Resubmission dreamer (NEW CRAN SUBMISSION)

This is a resubmission.  In this version, I have:

* Added \value to the .Rd files for the following functions:
** diagnostics.Rd
** dreamer_mcmc.Rd
** model_longitudinal.Rd
** model.Rd
** plot_trace.Rd
** post_perc_effect.Rd
** pr_med.Rd

* Resets par() values using on.exit() in
** plot_trace.dreamer
** plot_trace.dreamer_bma

## Test environments
* macOS-latest github continuous integration (devel and release)
* ubuntu-20.04 github continuous integration (devel and release)
* win-builder (devel and release)

## R CMD check results

There were no ERRORs or WARNINGs.

There was 1 NOTE from win-builder:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Richard Payne <paynestatistics@gmail.com>'
  
  New submission

## Downstream dependencies

There are currently no downstream dependencies for this package.
