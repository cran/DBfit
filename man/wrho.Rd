\name{wrho}
\alias{wrho}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating a new design matrix for Durbin stage 2
}
\description{
It returns a new design matrix for Durbin stage 2.
}
\usage{
wrho(xc, adjphi)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xc}{
centered design matrix, no column of ones
}
  \item{adjphi}{
(initial) estimate of rho in Durbin stage 1
}
}
\details{
see reference.
}
\references{
McKnight, S. D., McKean, J. W., and Huitema, B. E. (2000). A double bootstrap method
to analyze linear models with autoregressive error terms. Psychological methods,
5 (1), 87.
Shaofeng Zhang (2017). Ph.D. Dissertation.
}
\note{
This function is for internal use. The main function for users is \code{dbfit}.
}
