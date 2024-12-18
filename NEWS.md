# dreamer 3.2.0

* Remove automatic naming of models, now requires each model to be named (#34).
* Fix bug where candidate models could not be plotted with data (#38).

# dreamer 3.1.0

* Updated README.Rmd and vignettes (#1).
* Restore existing rjags modules on exit from `dreamer_mcmc()` (#4).
* Added tests to check `posterior()` with `predictive` argument (#6).
* Added "time" column for pr_eoi() output for longitudinal analyses (#9).
* Clean up jags seed and remove unused arguments on internal functions (#11).
* Update docs to make it clear sigma represents standard deviation (#14).
* Remove dreamer class -- use existing "dreamer_mcmc"" class instead (#17).
* Create safer class checks using `inherits()` (#16).
* Add missing examples to docs (#23).
* Adds additional data-type checks to `dreamer_mcmc()` (#3).
* Adds print() methods for models and mcmc outputs (#15).
* Uses temporary directories in test-plots.R.
* Update to README and vignette to show new print methods (#27).

# dreamer 3.0.0

* First open-source release.
