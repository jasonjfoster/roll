# roll <img src = "logo.png" align = "right" width = "120">

[![](https://api.travis-ci.org/jjf234/roll.svg)](https://travis-ci.org/jjf234/roll) [![](https://www.r-pkg.org/badges/version/roll)](https://cran.r-project.org/package=roll) [![](https://codecov.io/gh/jjf234/roll/graph/badge.svg)](https://codecov.io/github/jjf234/roll)
[![](https://cranlogs.r-pkg.org/badges/roll?color=brightgreen)](https://www.r-pkg.org/pkg/roll)

## Overview

`roll` is a package that provides fast and efficient computation of rolling and expanding statistics for time-series data.

The default algorithm in the `roll` package, and suitable for most applications, is an **online algorithm**. Based on the speed requirements and sequential nature of many problems in practice, online algorithms are a natural fit for computing rolling and expanding statistics of time-series data. That is, as observations are added and removed from a window, online algorithms update statistics and discard observations from memory; as a result, the amount of time to evaluate each function is significantly faster as the computation is independent of the window size (use the microbenchmark package to measure performance). In contrast, an offline algorithm requires all observations in memory to calculate the statistic for each window. Note that online algorithms are prone to loss of precision due to round-off error; hence, users can trade speed for accuracy and select the offline algorithm by setting the `online` argument to `FALSE`. Also, the RcppParallel package is used to parallelize the online algorithms across columns and across windows for the offline algorithms. 

As mentioned above, the numerical calculations use the RcppParallel package to parallelize rolling and expanding statistics of time-series data. The RcppParallel package provides a complete toolkit for creating safe, portable, high-performance parallel algorithms, built on top of the Intel Threading Building Blocks (TBB) and TinyThread libraries. By default, all the available cores on a machine are used for parallel algorithms. If users are either already taking advantage of parallelism or instead want to use a fixed number or proportion of threads, then set the number of threads in the RcppParallel package with the `RcppParallel::setThreadOptions` function.

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

Load the package and supply a dataset:

``` r
library(roll)

n <- 15
x <- rnorm(n)
y <- rnorm(n)
weights <- 0.9 ^ (n:1)
```

Then, to compute rolling and expanding means, use the `roll_mean` function:

``` r
# rolling means with complete windows
roll_mean(x, width = 5)

# rolling means with partial windows
roll_mean(x, width = 5, min_obs = 1)

# expanding means with partial windows
roll_mean(x, width = n, min_obs = 1)

# expanding means with weights and partial windows
roll_mean(x, width = n, weights = weights, min_obs = 1)
```

Also, the `roll_lm` function computes rolling and expanding regressions:

``` r
# rolling regressions with complete windows
roll_lm(x, y, width = 5)

# rolling regressions with partial windows
roll_lm(x, y, width = 5, min_obs = 1)

# expanding regressions with partial windows
roll_lm(x, y, width = n, min_obs = 1)

# expanding regressions with weights and partial windows 
roll_lm(x, y, width = n, weights = weights, min_obs = 1)
```

Note that the handling of missing values is supported as well (see the `min_obs`, `complete_obs`, and `na_restore` arguments).




