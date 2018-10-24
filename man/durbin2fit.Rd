\name{durbin2fit}
\alias{durbin2fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Durbin stage 2 fit
}
\description{
Function implements the Durbin stage 1 fit
}
\usage{
durbin2fit(yc, xc, adjphi, method, scores)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{yc}{
a transformed reponse variable
}
  \item{xc}{
a transformed design matrix
}
  \item{adjphi}{
the Durbin stage 1 estimate(s) of the autoregressive parameters rho
}
  \item{method}{
the method to be used for fitting. If "OLS", uses the ordinary least square; If "RANK", uses the rank-based fit.
}
 \item{scores}{
Default is Wilcoxon scores
}
}
\value{
\item{beta}{the estimates of regression coefficients}
\item{sigma}{the estimate of standard deviation of the white noise}
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


