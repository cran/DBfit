\name{durbin1xy}
\alias{durbin1xy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating New X and Y for Durbin Stage 1
}
\description{
Functions provides the tranformed reponse variable and model matrix for Durbin stage 1 fit. For details of the transformation, see the reference.
}
\usage{
durbin1xy(y, x, arp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
the orginal response variable
}
  \item{x}{
the orginal design matrix with first column of all one's (corresponding to the intercept)
}
  \item{arp}{
the order of autoregressive errors.
}
}
\references{
McKnight, S. D., McKean, J. W., and Huitema, B. E. (2000). A double bootstrap method
to analyze linear models with autoregressive error terms. Psychological methods,
5 (1), 87.
Shaofeng Zhang (2017). Ph.D. Dissertation.
}
