durbin1fit <-
function(y,x,arp,method,scores=scores){

     xy <- durbin1xy(y,x,arp)

     m <- ncol(xy)
     y2 <- xy[,m]
     x2 <- xy[,1:(m-1)]
     if(method=="OLS"){
       fit <- lm(y2 ~ x2)
     } else if (method=="RANK") {
       fit <- suppressWarnings(rfit(y2 ~ x2,scores=scores))
     }
     return(fit)
}
