\name{summary.dbfit}
\alias{summary.dbfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarize the double bootstrap (DB) fit
}
\description{
It summarizes the DB fit in a way that is similar to OLS \code{lm}.
}
\usage{
\method{summary}{dbfit}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a result of the call to \code{rfit}
}
  \item{\dots}{
additional arguments to be passed
}
}
\value{
\item{call}{the call to \code{rfit}}
\item{tab}{a table of point estimates, standard errors, t-ratios and p-values}
\item{rho1}{the Durbin two-stage estimate of rho}
\item{adjar}{the DB (final) estimate of rho}
\item{flag99}{an indicator; if 1, it indicates the original fit yields an estimate of rho to be 0.99. Only valid if order 1 is specified.
}
}
\examples{
# data(testdata)
# y<-testdata[,5]
# x<-testdata[,1:4]
# fit1<-dbfit(x,y,1) # OLS fit, default
# summary(fit1)
}