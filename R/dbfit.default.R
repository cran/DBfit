dbfit.default <-
function(x, y, arp, nbs=500, nbscov=500, conf=0.95, correction=TRUE,method="OLS",scores=Rfit::wscores, ...)
{
  x <- as.matrix(x)
  y <- as.numeric(y)
  est <- simula(x=x, y=y, arp=arp, nbs=nbs, nbscov=nbscov, conf=conf,method=method,scores=scores, ...)
  if ((arp == 1) && (est$adjar >= 0.99) && correction){
    est <- simulacorrection(x=x, y=y, arp=arp, nbs=nbs, nbscov=nbscov,method=method,scores=scores, ...)
  }
  est$call <- match.call()
  class(est) <- "dbfit"
  est
}
