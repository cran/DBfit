print.summary.dbfit <-
function(x, ...)
  {
    cat("Call:\n")
    print(x$call)
    
    cat("\nInitial rho:\n")
    print(x$rho1)
    
    cat("\nFinal rho:\n")
    print(c(x$adjar, x$CI_rho))
    
    cat("\nNonstationarity flag:\n")
    print(x$flag99)
    
    printCoefmat(x$tab, P.values = TRUE, has.Pvalue=TRUE)
  }
