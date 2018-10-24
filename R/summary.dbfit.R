summary.dbfit <-
function(object, ...)
{
  res <- list(call = object$call, tab = object$tabbeta, rho1 = unname(object$rho1), adjar = unname(object$adjar), 
               flag99 = object$flag99)
  class(res) <- "summary.dbfit"
  res
}
