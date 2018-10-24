\name{durbin1fit}
\alias{durbin1fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Durbin stage 1 fit
}
\description{
Function implements the Durbin stage 1 fit
}
\usage{
durbin1fit(y, x, arp, method, scores)
}
\arguments{
  \item{y}{
the response variable in stage 1, not the original response variable
}
  \item{x}{
the model matrix in stage 1, not the original design matrix
}
  \item{arp}{
the order of autoregressive errors.
}
  \item{method}{
the method to be used for fitting. If "OLS", uses the ordinary least square; If "RANK", uses the rank-based fit.
}
 \item{scores}{
Default is Wilcoxon scores
}
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

