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
Based on the speed requirements and sequential nature of many problems in practice, \bold{online algorithms} are a natural fit for computing rolling statistics of time-series data. That is, as observations are added and removed from a rolling window, online algorithms update statistics and discard observations from memory. The default algorithm in the \code{roll} package, and suitable for most applications, is an online algorithm; however, in some cases it is impossible to recover the information needed to update each statistic. Specifically, if the \code{weights} vector is an arbitrarily changing sequence then a standard algorithm is used instead to calculate the rolling statistic. In the former case, the algorithm is parallelized across columns via RcppParallel and across windows in the latter case. Note that online algorithms are prone to loss of precision due to round-off error; thus, in favor of accuracy over speed, the standard algorithm can be accessed by setting the \code{online} argument to \code{FALSE}.

As mentioned above, the numerical calculations use RcppParallel to parallelize rolling statistics of time-series data. RcppParallel provides a complete toolkit for creating safe, portable, high-performance parallel 
algorithms, built on top of the Intel Threading Building Blocks (TBB) and TinyThread libraries. By default, all the available cores on a machine are used for parallel algorithms. If users are 
either already taking advantage of parallelism or instead want to use a fixed number or proportion of 
threads, then set the number of threads in the RcppParallel package with the 
\code{\link[RcppParallel]{setThreadOptions}} function.
}
\author{
\packageAuthor{roll}
}
\references{
Welford, B.P. (1962). "Note on a method for calculating corrected sums of squares and products". \emph{Technometrics}, 4(3), 419–420.

West, D.H.D. (1979). "Updating Mean and Variance Estimates: An Improved Method". \emph{Communications of the ACM}, 22(9), 532-535.
}
