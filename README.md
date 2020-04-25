# roll <img src = "logo.png" align = "right" width = "120">

[![](https://api.travis-ci.org/jjf234/roll.svg)](https://travis-ci.org/jjf234/roll) [![](https://www.r-pkg.org/badges/version/roll)](https://cran.r-project.org/package=roll) [![](https://codecov.io/gh/jjf234/roll/graph/badge.svg)](https://codecov.io/github/jjf234/roll)
[![](https://cranlogs.r-pkg.org/badges/roll?color=brightgreen)](https://www.r-pkg.org/pkg/roll)

## Overview

`roll` is a package that provides fast and efficient computation of rolling and expanding statistics for time-series data.

The default algorithm in the `roll` package, and suitable for most applications, is an **online algorithm**. Based on the speed requirements and sequential nature of many problems in practice, online algorithms are a natural fit for computing rolling and expanding statistics of time-series data. That is, as observations are added and removed from a window, online algorithms update statistics and discard observations from memory; however, in some cases it is impossible to recover the information needed to update each statistic. Specifically, if the `weights` vector is an arbitrarily changing sequence then an offline algorithm is used instead to calculate the statistic. Also, in the former case, the algorithm is parallelized across columns via RcppParallel and across windows in the latter case. Note that online algorithms are prone to loss of precision due to round-off error; hence, users can trade speed for accuracy and select the offline algorithm by setting the `online` argument to `FALSE`.

As mentioned above, the numerical calculations use RcppParallel to parallelize rolling and expanding statistics of time-series data. RcppParallel provides a complete toolkit for creating safe, portable, high-performance parallel algorithms, built on top of the Intel Threading Building Blocks (TBB) and TinyThread libraries. By default, all the available cores on a machine are used for parallel algorithms. If users are either already taking advantage of parallelism or instead want to use a fixed number or proportion of threads, then set the number of threads in the RcppParallel package with the `RcppParallel::setThreadOptions` function.

## Installation

Get the released version from CRAN:

``` r
install.packages("roll")
```

Or the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("jjf234/roll")
```

## Usage

``` r
library(roll)

n <- 15
x <- rnorm(n)
y <- rnorm(n)
weights <- 0.9 ^ (n:1)

# -----------------------------------------------------------------------------

# rolling means with complete windows
roll_mean(x, width = 5)

# rolling means with partial windows
roll_mean(x, width = 5, min_obs = 1)

# expanding means with partial windows
roll_mean(x, width = n, min_obs = 1)

# expanding means with weights and partial windows
roll_mean(x, width = n, weights = weights, min_obs = 1)

# -----------------------------------------------------------------------------

# rolling regressions with complete windows
roll_lm(x, y, width = 5)

# rolling regressions with partial windows
roll_lm(x, y, width = 5, min_obs = 1)

# expanding regressions with partial windows
roll_lm(x, y, width = n, min_obs = 1)

# expanding regressions with weights and partial windows 
roll_lm(x, y, width = n, weights = weights, min_obs = 1)
```

Note that handling of missing values is also supported (see `min_obs`, `complete_obs`, and `na_restore` arguments).