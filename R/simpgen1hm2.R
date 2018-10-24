simpgen1hm2 <-
function (n1, n2, rho, beta = c(0, 0, 0, 0)) 
{
    n <- n1 + n2
    nstop <- 500 + n
    err <- rnorm(1)
    for (i in 2:nstop) {
        err[i] <- rho * err[i - 1] + rnorm(1)
    }
    errs <- err[501:nstop]
    xmat <- hmdesign2(n1, n2)
    y <- xmat %*% beta + errs
    mat <- cbind(xmat, y)
    return(mat)
}
