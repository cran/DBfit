\name{hmdesign2}
\alias{hmdesign2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
the Two-Phase Design Matrix
}
\description{
Returns the design matrix for a two-phase intervention model.
}
\usage{
hmdesign2(n1, n2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n1}{
number of obs in phase 1
}
  \item{n2}{
number of obs in phase 2
}
}
\details{
It returns a matrix of 4 columns.
As discussed in Huitema, Mckean, & Mcknight (1999), in
two-phase design: beta0 = intercept, beta1 = slope for Phase 1, beta2 = level change from Phase
1 to Phase 2, and beta3 slope change from Phase 1 to Phase 2.
}

\references{
Huitema, B. E., Mckean, J. W., & Mcknight, S. (1999). Autocorrelation effects on least-
squares intervention analysis of short time series. Educational and Psychological
Measurement, 59 (5), 767-786.
}

\examples{
n1 <- 15
n2 <- 15
hmdesign2(n1, n2)
}