# roll <img src = "man/figures/logo.png" align = "right" width = "120">

[![](https://github.com/jasonjfoster/roll/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/jasonjfoster/roll/actions/workflows/check-standard.yaml)
[![](https://www.r-pkg.org/badges/version/roll)](https://cran.r-project.org/package=roll)
[![](https://codecov.io/gh/jasonjfoster/roll/graph/badge.svg)](https://app.codecov.io/github/jasonjfoster/roll)
[![](https://cranlogs.r-pkg.org/badges/roll?color=brightgreen)](https://www.r-pkg.org/pkg/roll)

## Overview

'roll' provides fast and efficient computation of rolling and expanding statistics for time-series data.

The default algorithm in the 'roll' package is an online algorithm that, as observations are added to and removed from a window, updates statistics and discards observations from memory (Welford, 1962, <doi:10.1080/00401706.1962.10490022>; West, 1979, <doi:10.1145/359146.359153>); as a result, the amount of time to evaluate each function is significantly shorter as the computation is independent of the window. In contrast, an offline algorithm requires all observations in memory to calculate the statistic for each window, so users can trade speed for accuracy and select the offline algorithm by setting the `online` argument to `FALSE`. Quantiles are computed from the inverse of the empirical distribution function with averaging at discontinuities (Hyndman and Fan, 1996, <doi:10.1080/00031305.1996.10473566>). Use cases include:

* **Rolling summary statistics**: means, standard deviations, medians, and quantiles that adapt to the most recent window of observations
* **Time-varying relationships**: rolling variances, covariances, correlations, and linear regressions between variables
* **Feature engineering**: rolling z-scores and standardized series for forecasting and signal construction

The package supports rolling and expanding windows, weights, and handling of missing values via the min_obs, complete_obs, and na_restore arguments. The implementation uses 'RcppParallel' to parallelize the online algorithms across columns and the offline algorithms across windows.

## Installation

Install the released version from CRAN:

```r
install.packages("roll")
```

Or the development version from GitHub:

```r
# install.packages("pak")
pak::pak("jasonjfoster/roll")
```

## Usage

Load the package and supply a dataset:

```r
library(roll)

n <- 15
x <- rnorm(n)
y <- rnorm(n)
weights <- 0.9 ^ (n:1)
```

Then, to compute rolling and expanding means, use the `roll_mean()` function:

```r
# rolling means with complete windows
roll_mean(x, width = 5)

# rolling means with partial windows
roll_mean(x, width = 5, min_obs = 1)

# expanding means with partial windows
roll_mean(x, width = n, min_obs = 1)

# expanding means with partial windows and weights
roll_mean(x, width = n, min_obs = 1, weights = weights)
```

Or use the `roll_lm()` function to compute rolling and expanding regressions:

```r
# rolling regressions with complete windows
roll_lm(x, y, width = 5)

# rolling regressions with partial windows
roll_lm(x, y, width = 5, min_obs = 1)

# expanding regressions with partial windows
roll_lm(x, y, width = n, min_obs = 1)

# expanding regressions with partial windows and weights
roll_lm(x, y, width = n, min_obs = 1, weights = weights)
```

Handling of missing values is also supported (see the `min_obs`, `complete_obs`, and `na_restore` arguments).

## References

Hyndman, R.J. and Fan, Y. (1996). "Sample Quantiles in Statistical Packages." *The American Statistician* 50 (4): 361-365. <doi:10.1080/00031305.1996.10473566>

Welford, B.P. (1962). "Note on a Method for Calculating Corrected Sums of Squares and Products." *Technometrics* 4 (3): 419-420. <doi:10.1080/00401706.1962.10490022>

West, D.H.D. (1979). "Updating Mean and Variance Estimates: An Improved Method." *Communications of the ACM* 22 (9): 532-535. <doi:10.1145/359146.359153>
