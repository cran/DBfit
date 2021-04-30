\name{hypothmat}
\alias{hypothmat}
\title{
General Linear Tests of the regression coefficients
}
\description{
Performs general linear tests of the regressio coefficients.
}
\usage{
hypothmat(sfit, mmat, n, p)
}
\arguments{
  \item{sfit}{
the result of a call to dbfit.
}
  \item{mmat}{
a full row rank q*(p+1) matrix, where q is the row number of the matrix and p is number of independent variables.
}
  \item{n}{
total number of observations.
}
  \item{p}{
number of independent variables.
}
}
\details{
This functions performs the general linear F-test of the form H0: Mb = 0 vs HA: Mb != 0.
}
\value{
\item{tst}{the test statistic}
\item{pvf}{the p-value of the F-test}
}
\references{
McKnight, S. D., McKean, J. W., and Huitema, B. E. (2000). A double bootstrap method
to analyze linear models with autoregressive error terms. Psychological methods,
5 (1), 87.
Shaofeng Zhang (2017). Ph.D. Dissertation.
}
\examples{
# data(testdata)
# y<-testdata[,5]
# x<-testdata[,1:4]
# fit1<-dbfit(x,y,1) # OLS fit, default
# a test that H0: b1 = b3 vs HA: b1 != b3
# mat<-matrix(c(1,0,0,-1),nrow=1) 
# hypothmat(sfit=fit1,mmat=mat,n=40,p=4)
}

