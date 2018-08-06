## News for package 'roll'

### Changes in roll version 1.1.0 (YYYY-MM-DD)

* Added `online` argument to process observations using an online algorithm by default

* `roll_lm` function now returns standard errors ([#7](https://github.com/jjf234/roll/issues/7))

* Simplified checks for `width` and `min_obs` arguments ([#3](https://github.com/jjf234/roll/issues/3))

* Added `y` argument to `roll_cov` and `roll_cor` functions ([#2](https://github.com/jjf234/roll/issues/2))

* Updated `src/Makevars` and `src/Makevars.win` files to what the `RcppArmadillo` skeleton default now uses to more fully utilize OpenMP

    * Note: if users take advantage of parallelism using multithreaded libraries, then limit the number of cores in the `RcppParallel` package to one with the `setThreadOptions` function
    
* Deprecated less common functions (`roll_eigen`, `roll_vif`, `roll_pcr`) and arguments (`scale` and `center` in the `roll_lm` function); also removed the `parallel_for` argument in favor of a new approach used internally

### Changes in roll version 1.0.7 (2017-05-01)

* New `roll_sum` and `roll_prod` functions for computing rolling sums and products, respectively, of time-series data

* Added `init.c` file with calls to `R_registerRoutines()` and `R_useDynamicSymbols()`; also uses `.registration = TRUE` in `useDynLib` in `NAMESPACE`

### Changes in roll version 1.0.6 (2016-09-19)

* Added `intercept` argument to `roll_lm` and `roll_pcr` functions

* Turned on `CXXSTD = CXX11` to enforce adherence to the C++11 standard

### Changes in roll version 1.0.5 (2016-08-29)

* Added a section on examples to the `README` file

* Fixed an issue in the `src/Makevars` and `src/Makevars.win` files ([#1](https://github.com/jjf234/roll/issues/1))

### Changes in roll version 1.0.4 (2016-07-05)

* `roll_lm` and `roll_pcr` functions have been enhanced:

    * `y` can now be a matrix or xts object with multiple dependent variables

    * Added shorthand arguments for `center` and `scale`

* New `roll_scale` function for computing rolling scaling and centering of time-series data