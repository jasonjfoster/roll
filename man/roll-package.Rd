\name{roll-package}
\alias{roll-package}
\alias{roll}
\docType{package}
\title{
\packageTitle{roll}
}
\description{
\packageDescription{roll}
}
\details{
\code{roll} is a package that provides fast and efficient computation of rolling and expanding statistics for time-series data.

The default algorithm in the \code{roll} package, and suitable for most applications, is an \bold{online algorithm}. Based on the speed requirements and sequential nature of many problems in practice, online algorithms are a natural fit for computing rolling and expanding statistics of time-series data. That is, as observations are added and removed from a window, online algorithms update statistics and discard observations from memory; however, in some cases it is impossible to recover the information needed to update each statistic. Specifically, if the \code{weights} vector is an arbitrarily changing sequence then an offline algorithm is used instead to calculate the statistic. For reference, an offline algorithm requires all observations in memory to calculate the statistic for each window. Note that online algorithms are prone to loss of precision due to round-off error; hence, users can trade speed for accuracy and select the offline algorithm by setting the \code{online} argument to \code{FALSE}. Also, the online algorithm is parallelized across columns via RcppParallel and across windows for offline algorithms. 

As mentioned above, the numerical calculations use RcppParallel to parallelize rolling and expanding statistics of time-series data. RcppParallel provides a complete toolkit for creating safe, portable, high-performance parallel  algorithms, built on top of the Intel Threading Building Blocks (TBB) and TinyThread libraries. By default, all the available cores on a machine are used for parallel algorithms. If users are  either already taking advantage of parallelism or instead want to use a fixed number or proportion of threads, then set the number of threads in the RcppParallel package with the \code{RcppParallel::setThreadOptions} function.
}
\author{
\packageAuthor{roll}
}
\references{
Welford, B.P. (1962). "Note on a Method for Calculating Corrected Sums of Squares and Products." \emph{Technometrics}, 4(3), 419-420.

West, D.H.D. (1979). "Updating Mean and Variance Estimates: An Improved Method." \emph{Communications of the ACM}, 22(9), 532-535.
}
