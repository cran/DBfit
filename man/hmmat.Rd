\name{hmmat}
\alias{hmmat}
\title{
K-Phase Design Matrix
}
\description{
Returns the design matrix for a general k-phase intervention model
}
\usage{
hmmat(vecss, k)
}
\arguments{
  \item{vecss}{
a vector of length k with each element being the number of observations in each phase
}
  \item{k}{
number of phases
}
}
\details{
It returns a matrix of 2*k columns. The design can be unbalanced, i.e. each phase has different observations.
}
\references{
Huitema, B. E., Mckean, J. W., & Mcknight, S. (1999). Autocorrelation effects on least-
squares intervention analysis of short time series. Educational and Psychological
Measurement, 59 (5), 767-786.
}
\seealso{
\code{\link{hmdesign2}}
}
\examples{
# a three-phase design matrix
hmmat(c(10,10,10),3)
}
