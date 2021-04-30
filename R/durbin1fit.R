durbin1fit <-
function(y,x,arp,method,scores=scores){

     xy <- durbin1xy(y,x,arp)

     m <- ncol(xy)
     y2 <- xy[,m]
     x2 <- xy[,1:(m-1)]
     if(method=="OLS"){
       fit <- lm(y2 ~ x2)
     } else if (method=="RANK") {
	     p2 <- ncol(x) - 1
	     n2 <- nrow(xy)
	     x3 <- fullr(cbind(rep(1,n2),x2),1+arp+p2)
#       fit <- suppressWarnings(rfit(y2 ~ x2,scores=scores))
        fit <- rfit(y2 ~ x3[,-1],scores=scores)
     }
     return(fit)
}
