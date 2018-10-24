\name{lagx}
\alias{lagx}
\alias{lagmat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Lag Functions
}
\description{
For preparing the transformed x and y in the Durbin stage 1 fit
}
\usage{
lagx(x, s1, s2)
lagmat(x, p)
}
\arguments{
  \item{x}{
a vector or the design matrix, including intercept, i.e. the first column being ones.
}
  \item{s1}{
starting index of the slice.
}
  \item{s2}{
end index of the slice.
}
  \item{p}{
the order of autoregressive errors.
}
}
\note{
These function are for internal use.
}

