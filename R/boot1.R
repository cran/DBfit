boot1 <-
function(y,phi1,arp,nbs,x,allb,method,scores=scores){

     upper <- .99
     lower <- -.99
     for(j in 1:arp){
         if(phi1[j] < lower){phi1[j] <- lower}
         if(phi1[j] > upper){phi1[j] <- upper}
     }

     icent <- 0
     phia <- phi1
     xcopy <- x
     xcpy <- x
     n <- length(y)
     yy <- y
     p <-length(x[1,])
     adj <- (n - arp - p)/(n - arp - 2*p)

     ypart <- nurho(y,phi1)
     n1 <- n - arp

     xpart <- wrho(x,phi1)
     ehat <- ypart - xpart%*%allb

     ehat <- (ehat - mean(ehat))*sqrt(adj)

#     bootstrap loop

     ind <- 1:n1
     bsr1 <- matrix(rep(0,nbs*arp),ncol=arp)

     for(nbk in 1:nbs){
         ystar <- yy[1:arp]
         ind2 <- sample(ind,n1,replace=TRUE)
         estar <- ehat[ind2]

         for(i in (arp+1):n){
              ypart <- 0
              for(k in 1:arp){
                  ypart <- ypart + phia[k]*ystar[i-k]
              }
              xpart <- 0
              for(j in 1:p){
                  for(k in 1:arp){
                       if(k == 1){
                            xpart <- xpart + allb[j]*(xcopy[i,j]-phia[k]*xcopy[i-k,j])
                       } else {
                            xpart <- xpart - allb[j]*phia[k]*xcopy[i-k,j]
                       }
                  }
              }
              ystar[i] <- ypart + xpart + estar[i-arp]
         }
         if(icent == 1){
              avey <- mean(ystar)
              ystar <- ystar - avey
         }
         dfit1 <- durbin1fit(ystar,xcpy,arp,method=method,scores=scores)
         adjphi1 <- dfit1$coef[2:(arp+1)]
         for(k in 1:arp){
               if(adjphi1[k] < lower){adjphi1[k] <- lower}
               if(adjphi1[k] > upper){adjphi1[k] <- upper}
               bsr1[nbk,k] <- phi1[k] - adjphi1[k]
         }
     }

     bias1 <- apply(bsr1,2,mean)
     return(bias1)
}
