\name{boot2}
\alias{boot2}
\title{
First Boostrap Procedure For parameter estimations
}
\description{
Function performing the second bootstrap procedure to yield the inference of the regression coefficients
}
\usage{
boot2(y, xcopy, phi1, beta, nbs, method, scores)
}
\arguments{
  \item{y}{
the response variable
}
  \item{xcopy}{
the original design matrix (including intercept), without centering
}
  \item{phi1}{
the estimate of the autoregressive parameter rho from the first bootstrap procedure
}
  \item{beta}{
the estimates of the regression coefficients from the first bootstrap procedure
}
  \item{nbs}{
the bootstrap size
}
  \item{method}{
If "OLS", uses the ordinary least square; If "RANK", uses rank-based fit
}
 \item{scores}{
Default is Wilcoxon scores
}
}
\value{
\item{betacov}{
the estimate of var-cov matrix of betas
}
\item{allbeta}{
the estimates of betas inside of the second bootstrap, not the final estimates of betas. The final estimates of betas are still from \code{boot1}.
}
\item{rhostar}{
the estimates of rho inside of the second bootstrap, not the final estimates of rho. The final estimate(s) of rho are still from \code{boot1}.
}
\item{MSEstar}{
MSE used inside of the second bootstrap.
}
}
\note{
This function is for internal use. The main function for users is \code{dbfit}
}



