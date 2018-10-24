\name{boot1}
\alias{boot1}
\title{
First Boostrap Procedure For parameter estimations
}
\description{
Function performing the first bootstrap procedure to yield the parameter estimates
}
\usage{
boot1(y, phi1, arp, nbs, x, allb, method, scores)
}
\arguments{
  \item{y}{
the response variable
}
  \item{phi1}{
the Durbin two-stage estimate of the autoregressive parameter rho
}
  \item{arp}{
the order of autoregressive errors
}
  \item{nbs}{
the bootstrap size
}
  \item{x}{
the original design matrix (including intercept), without centering
}
  \item{allb}{
all the Durbin two-stage estimates of the regression coefficients 
}
  \item{method}{
If "OLS", uses the ordinary least square; If "RANK", uses the rank-based fit
}
 \item{scores}{
Default is Wilcoxon scores
}
}
\value{
An estimate of the bias is returned
}
\note{
This function is for internal use. The main function for users is \code{dbfit}.
}
