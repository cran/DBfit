\name{simpgen1hm2}
\alias{simpgen1hm2}
\title{
Simulation Data Generating Function
}
\description{
Generates the simulation data for a two-phase intervention model.
}
\usage{
simpgen1hm2(n1, n2, rho, beta = c(0, 0, 0, 0))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n1}{
number of obs in phase 1
}
  \item{n2}{
number of obs in phase 2
}
  \item{rho}{
pre-defined autoregressive parameter(s)
}
  \item{beta}{
pre-defined regression coefficients
}
}
\details{
This function is used for simulations when developing the package. With pre-defined sample sizes in both phases and parameters, it returns a simulated data.
}
\value{
\item{mat}{a matrix containing the simulation data. The last column is the response variable. All other columns make up the design matrix.}
}
\seealso{
\code{\link{hmdesign2}}
}
\examples{
 n1 <- 15
 n2 <- 15
 rho <- 0.6
 beta <- c(0,0,0,0)
 dat <- simpgen1hm2(n1, n2, rho, beta)
 dat
}

