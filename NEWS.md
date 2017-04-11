## News for package 'roll'

### Changes in roll version 1.0.7 (YYYY-MM-DD)

* New `roll_sum` and `roll_prod` functions for computing rolling sums and products, respectively, of time-series data

* Added `init.c` files with calls to `R_registerRoutines()` and `R_useDynamicSymbols()`; also uses `.registration = TRUE` in `useDynLib` in `NAMESPACE`

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