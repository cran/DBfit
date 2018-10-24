dbfit.formula <-
function(formula, data=list(),arp,nbs=500,nbscov=500,conf=0.95,correction=TRUE,method="OLS",scores=scores,...)
{
  mf <- model.frame(formula=formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  est <- dbfit.default(x=x, y=y,arp=arp,nbs=nbs,nbscov=nbscov,conf=conf,correction=correction,method=method,scores=scores, ...)
  est$call <- match.call()
  est$formula <- formula
  est
}
