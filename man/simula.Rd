\name{simula}
\alias{simula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Work Horse Function to implement the Double Bootstrap method
}
\description{
\code{simula} is the original work horse function to implement the DB method. However, when this function returns an estimate of rho to be .99, another work horse function \code{simulacorrection} kicks in.
}
\usage{
simula(x, y, arp, nbs, nbscov, conf, method, scores)
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
  \item{conf}{
the confidence level of CI for rho, default is 0.95. 
}
  \item{method}{
the method to be used for fitting. If "OLS", uses the ordinary least square \code{lm}; If "RANK", uses the rank-based fit \code{rfit}.
}
 \item{scores}{
Default is Wilcoxon scores
}
}
\details{
see \code{\link{dbfit}}.
}
\references{
McKnight, S. D., McKean, J. W., and Huitema, B. E. (2000). A double bootstrap method
to analyze linear models with autoregressive error terms. Psychological methods,
5 (1), 87.
Shaofeng Zhang (2017). Ph.D. Dissertation.
}
\note{
Users should use \code{dbfit} to perform the analysis.
}
\seealso{
\code{\link{dbfit}}.
}

