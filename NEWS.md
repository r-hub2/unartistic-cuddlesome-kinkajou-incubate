# incubate 1.2.1
* fix package help: due to changes in roxygen

# incubate 1.2.0
* check for minimal number of observations and fail early if not enough observations
* `test_diff`:
    * test-statistic gets lower bound of 0 enforced (a restricted model can not have better fit than unrestricted model)
    * deactivate Anderson-Darling (AD) GOF-test as its performance in simulations under two-group setting was unsatisfactory
* include simulation R-scripts in package under `inst/scripts/` folder. When the package is installed the scripts are found at in the -directory `scripts/` within the package.
* rename methods:
    * 'MSE' => 'MPSE' (as the term MSE has already other meanings)
    * 'MLE' => 'MLE0' (to indicate that this is the standard MLE that is not appropriate for delay models)
* rework implementation of objective-function: do not rely on attributes but use the environment of the objective function instead. This means fewer copies of information are made.

# incubate 1.1.9
* New implementation of delay smoothing: use objective function along candidate delay values to guide in smoothing
    * define an interval for delay candidate values
    * sample from this interval according to the objective function: the better the value of the objective function the more likely it is to sample a value. For speed of computation, we vary only delay and keep the other parameters at their estimate.
    * use `smd_factor` to either concentrate or flatten around the best delay estimate
* Better implementation for MLE-fitting:
    * Support MLE-fitting also for Weibull (one group setting only)
    * Move MLE-objective function into a general objective function factory method

# incubate 1.1.8
* correct URL for bug reports in DESCRIPTION

# incubate 1.1.7
* more robust test that passes also on Apple-silicon
* doc: add homepage in DESCRIPTION

# incubate 1.1.6
* use LGPL-licence

# incubate 1.1.5
* doc: further improvements (document all return values)

# incubate 1.1.4
* doc: small improvements

# incubate 1.1.3
* README: better intro text and more text for example

# incubate 1.1.2
* small changes for CRAN-warnings
* test: bigger data set for delay estimation

# incubate 1.1.1
* Fix errors/warnings from check
    * document every function parameter
    * S3-methods have right parameters
    * specify package name for external functions (e.g., from `stats`)
    * global variables (use `.data` and `aes_`)
* add `cran-comments.md` file

# incubate 1.1.0
* confint: use better default values
* confint: when using bootstrap smoothing for delay, make sure the smoothing keeps limited for tiny estimated rates lambda (see `SMD_MINRATE`)
* bug fix: `plot.incubate_fit` has legible colour legend title

# incubate 1.0.0
* data: add experimental data `stankovic` from a publication Stankovic et al EGFL7 enhances surface expression of integrin doi:10.15252/emmm.201708420 (2018) as example data set
* This version of package is also published at <https://zenodo.org/record/6462382>

# incubate 0.9.9
* `power_diff`: improve heuristics in iterative search for sample size, given power

# incubate 0.9.8
* `test_diff`: change name of P-values (within object)
* `power_diff`: allow to specify power and do iterative computation

# incubate 0.9.7
* confint: smooth delay for parametric bootstrap & delayed exponential model via rectified Gaussian (instead of normal Gaussian)

# incubate 0.9.6
* confint: add option to smooth delay for parametric bootstrap & delayed exponential model

# incubate 0.9.5
* confint: redefine `logshift` to mean what is *added* to the argument of log, before it was subtracted. This meaning seems more natural.

# incubate 0.9.4
* confint: logshift accepts also `NA` and `NULL`. Useful for simulation studies with optional logshift-setting.

# incubate 0.9.3
* confint: parameterized logshift for delay, new default is 0.01

# incubate 0.9.2
* confint: log-transformation revisited: shift mainly for delay (and hardly a shift for the other parameters)
* confint: extend test coverage to confint routines

# incubate 0.9.1
* confint: allow for log-transformation of bootstrapped coefficients

# incubate 0.9
* confint: separate bootstrap data generation and inference step. This allows for more efficient confint-simulations.

# incubate 0.8
* enhance test-coverage for package
* reorganize tests to better match the R-script

# incubate 0.7.6
* `test_diff`: do GOF-tests also for unrestricted model (e.g. `gof_pearson1`) besides for the restricted null model (renamed to `gof_pearson0`)

# incubate 0.7.5
* remove dependency on `dplyr`

# incubate 0.7.4
* add `bs_infer = 't0'` 
* add +3 to degrees of freedom for t-quantiles in `bs_infer='t'` and `='t0'`
* version bump belated: these features were temporarily also released as 0.7.3

# incubate 0.7.3
* fixed mapping of bootstrap inference names for `boot`-package

# incubate 0.7.2
* allow to use `boot` package to calculate confidence interval

# incubate 0.7.1
* allow to choose which tests to perform in `test_diff` (helpful when calculation of AD GOF-tests fails on older R-installations)

# incubate 0.7
* change of S3-classes:
    * for delay model fitting the new class is `incubate_fit`
    * for tests the new class is `incubate_test`
* confidence intervals supports now simple data generation to draw from the data with replacement (besides parametric bootstrap)
* preliminary support for MLE-fitting: this allows to compare confidence intervals based on MSE vs MLE

# incubate 0.6.1
* make handling of ties in calculation of spacings within the objective function more robust

# incubate 0.6
* implement AD-GOF test within `test_GOF`

# incubate 0.5
* add separate function `test_GOF` for goodness of fit (GOF) tests based on a fitted model
* add Moran's GOF-test

# incubate 0.4
* bug fix in MSE-criterion
* bug fix in tie-handling
* function `estimRoundingError` to estimate rounding error (for tie handling)

# incubate 0.3
* simplified power function (dropped stuff and conventions from `sscn`-package)

# incubate 0.2
* `ties='equidist'` is default method to handle ties
  
# incubate 0.1
* more choice in handling of ties: `equidist` and `random`
* rename function `test_delay_diff` to `test_diff` 

# incubate 0.0.5
* add confidence intervals (based on basic bootstrap) for model parameters

# incubate 0.0.4
* Weibull: `bind=` is implemented for all possible parameter combinations
* internal clean up

# incubate 0.0.3
* more warnings when MSE-optimization fails
* bug fix for power simulation
* bug fix in testing

# incubate 0.0.2
* initial start as a move out from package `sscn` (git history is not moved along, though)
* code robustness: add `try` code at lower level, directly at where `optim` is called.

