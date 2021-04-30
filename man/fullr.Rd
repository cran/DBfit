\name{fullr}
\alias{fullr}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 QR decomposition for non-full rank design matrix for Rfit.
}
\description{
With Rfit recent update, it cannot return partial results with sigular design matrix (as opposed to lm). This function uses QR decomposition for Rfit to resolve this issue, so that dbfit can run robust version.
}
\usage{
fullr(x, p1)

}
\arguments{
  \item{x}{
design matrix, including intercept, i.e. the first column being ones.
}
  \item{p1}{
number of first few columns of x that are lineraly independent.
}
}
\note{
This function is for internal use.
}

