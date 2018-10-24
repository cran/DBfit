\name{rhoci2}
\alias{rhoci2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A fisher type CI of the autoregressive parameter rho
}
\description{
This function returns a Fisher type CI for rho, which is then used to correct the .99 cases.
}
\usage{
rhoci2(n, rho, cv)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
total number of observations
}
  \item{rho}{
final estimate of rho, usually .99.
}
  \item{cv}{
critical value for CI
}
}
\details{
see reference.
}
\references{
Shaofeng Zhang (2017). Ph.D. Dissertation.
Rao, C. R. (1952). Advanced statistical methods in biometric research. p. 231
}
\note{
This function is for internal use.
}

