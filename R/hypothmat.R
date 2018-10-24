hypothmat <-
function(sfit,mmat,n,p){

      q <- length(mmat[,1])
      bpart <- mmat%*%sfit$coefficients
      varpart <- mmat%*%sfit$betacov%*%t(mmat)
      tst <- t(bpart)%*%solve(varpart)%*%bpart
      # pv <- 1 - pchisq(tst,q)
      pvf <- 1 - pf(tst/q,q,n-p)
      hypothmat <- c(tst,pvf)
      return(hypothmat)
}
