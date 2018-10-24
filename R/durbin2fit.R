durbin2fit <-
function(yc,xc,adjphi,method,scores=scores){
     p<-ncol(xc)
     arp <- length(adjphi)
     n <- length(yc)

     nuy <- nurho(yc,adjphi)
     wx <- wrho(xc,adjphi)
     
     if(method=="OLS"){
       fitls <- lm(nuy ~ wx - 1)
     } else if (method=="RANK") {
       fitls <- rfit(nuy ~ wx - 1,scores=scores)
     }
     
     resd <- fitls$resid
     beta <- fitls$coef
     n1 <- length(resd)
     
     sigma2 <- var(resd)
     sigma2 <- (n1/(n1 - 2*p -1))*sigma2
     sigma <- sqrt(sigma2)

     list(beta=beta,sigma=sigma)
}
