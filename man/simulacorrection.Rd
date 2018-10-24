\name{simulacorrection}
\alias{simulacorrection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Work Horse Function to Implement the Double Bootstrap Method For .99 Cases 
}
\description{
When function \code{simula} returns an estimate of rho to be .99, this function kicks in and ouputs a corrected estimate of rho. Currently, this only works for order 1, i.e. for order > 1, this correction will not get involved.
}
\usage{
simulacorrection(x, y, arp, nbs, nbscov, method, scores)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
the design matrix, including intercept, i.e. the first column being ones.
}
  \item{y}{
the response variable.
}
  \item{arp}{
the order of autoregressive errors.
}
  \item{nbs}{
the bootstrap size for the first bootstrap procedure. Default is 500.
}
  \item{nbscov}{
the bootstrap size for the second bootstrap procedure. Default is 500.
}
  \item{method}{
the method to be used for fitting. If "OLS", uses the ordinary least square \code{lm}; If "RANK", uses the rank-based fit \code{rfit}.
}
 \item{scores}{
Default is Wilcoxon scores
}
}
\details{
If 0.99 problem is detected, then construct Fisher
CI for both initial estimate (in Durbin stage 1) and first bias-corrected estimate (perform only one bootstrap, instead of a loop); if the midpoint of latter is smaller than 0.95, then this midpoint is the final
estimate for rho; otherwise the midpoint of the former CI is the final estimate.

By default, when function \code{simula} returns an estimate of rho to be .99, this function kicks in and ouputs a corrected estimate of rho. However, users can turn the auto correction off by setting correction="FALSE" in \code{dbfit}. Users are encouraged to investigate why the stationarity assumption is violated based on their experience of time series analysis and knowledge of the data.
}

\references{
Shaofeng Zhang (2017). Ph.D. Dissertation.
}
\note{
Users should use \code{dbfit} to perform the analysis.
}
\seealso{
\code{\link{dbfit}}.
}