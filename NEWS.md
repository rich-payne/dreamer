# dreamer 3.0.0.9000

* Updated README.Rmd and vignettes (#1).
* Restore existing rjags modules on exit from `dreamer_mcmc()` (#4).
* Added tests to check `posterior()` with `predictive` argument (#6).
* Added "time" column for pr_eoi() output for longitudinal analyses (#9).
* Clean up jags seed and remove unused arguments on internal functions (#11).
* Update docs to make it clear sigma represents standard deviation (#14).
* Remove dreamer class -- use existing "dreamer_mcmc"" class instead (#17).
* Create safer class checks using `inherits()` (#16).

# dreamer 3.0.0

* First open-source release.
