# roll

## Version 1.1.5

* Added hex sticker and improved documentation ([#27](https://github.com/jjf234/roll/issues/27), [#28](https://github.com/jjf234/roll/issues/28))

* New `roll_quantile` function for computing rolling and expanding quantiles of time-series data

    * Note: `roll_quantile` function is not calculated using an online algorithm

## Version 1.1.4

* Fixed issues related to floating point arithmetic ([#23](https://github.com/jjf234/roll/issues/23), [#24](https://github.com/jjf234/roll/issues/24))

## Version 1.1.3

* Added vector support to each function ([#20](https://github.com/jjf234/roll/issues/20))

* Implemented efficient algorithms for `roll_min` and `roll_max` functions

* New `roll_idxmin` and `roll_idxmax` functions for computing rolling and expanding indices of minimums and maximums, respectively, of time-series data ([#22](https://github.com/jjf234/roll/issues/22))

## Version 1.1.2

* New `roll_median`, `roll_min`, `roll_max`, `roll_any`, and `roll_all` functions for computing rolling and expanding medians, minimums, maximums, any, and all, respectively, of time-series data ([#4](https://github.com/jjf234/roll/issues/4), [#13](https://github.com/jjf234/roll/issues/13), [#14](https://github.com/jjf234/roll/issues/14))
    
    * Note: `roll_median`, `roll_min`, and `roll_max` functions are not calculated using online algorithms

## Version 1.1.1

* Added `online` argument to process observations using online algorithms by default

* `roll_lm` function now returns standard errors ([#7](https://github.com/jjf234/roll/issues/7))

* Simplified checks for `width` and `min_obs` arguments ([#3](https://github.com/jjf234/roll/issues/3))

* Added `y` argument to `roll_cov` and `roll_cor` functions ([#2](https://github.com/jjf234/roll/issues/2))

* Deprecated less common functions (`roll_eigen`, `roll_vif`, and `roll_pcr`) and arguments (`scale` and `center` in the `roll_lm` function); also removed the `parallel_for` argument in favor of a new approach used internally

## Version 1.0.7

* New `roll_sum` and `roll_prod` functions for computing rolling and expanding sums and products, respectively, of time-series data

## Version 1.0.6

* Added `intercept` argument to `roll_lm` and `roll_pcr` functions

## Version 1.0.5

* Added section on examples to the `README` file

* Fixed an issue in the `src/Makevars` and `src/Makevars.win` files ([#1](https://github.com/jjf234/roll/issues/1))

## Version 1.0.4

* `roll_lm` and `roll_pcr` functions have been enhanced:

    * `y` can now be a matrix or xts object with multiple dependent variables

    * Added shorthand arguments for `center` and `scale`

* New `roll_scale` function for computing rolling and expanding scaling and centering of time-series data