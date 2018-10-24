rhoci2 <-
function(n,rho,cv){

    rat <- (1-rho)/(1+rho)
    u1 <- exp(-cv*(2/sqrt(n-3)))
    u2 <- exp(cv*(2/sqrt(n-3)))
    ub <- (1-rat*u1)/(1+rat*u1)
    lb <- (1-rat*u2)/(1+rat*u2)
    rhoci <- c(rho,lb,ub)
    return(rhoci)
}
