boot2 <-
function(y,xcopy,phi1,beta,nbs,method,scores=scores){
#
#   beta includes intercept

      arp <- length(phi1)
      phi <- phi1
      p <- length(beta)
      n <- length(y)
      zn1 <- n-p-arp
      zn2 <- n-2*p-arp
      if(zn1 < 0){zn1 <- 1}
      if(zn2 < 0){zn2 <- 1}
      adj <- zn1/zn2
      icent <- 1

      n1 <- n - arp
      ones <- rep(1,n)
      proj1 <- ones%*%t(ones)/n
      x2 <- as.matrix(xcopy[,2:p])
      xbar <- apply(x2,2,mean)
      xc <- xcopy[,2:p] - proj1%*%xcopy[,2:p]
      x <- xc
      p <- p - 1
      ehat <- rep(0,n1)

      for(i in (arp+1):n){
          ypart <- y[i]
          for(k in 1:arp){ypart <- ypart - phi[k]*y[i-k]}
          xpart <- 0
          for(j in 1:(p+icent)){
               for(k in 1:arp){
                   if(k == 1){
                        xpart <- xpart + beta[j]*(xcopy[i,j] - phi[k]*xcopy[i-k,j])
                   } else {
                       xpart <- xpart - beta[j]*phi[k]*xcopy[i-k,j]
                   }
               }
          }
          ehat[i-arp] <- ypart - xpart
       }
       ehat <- (ehat - mean(ehat))*sqrt(adj)
       sigehat <- sd(ehat)*sqrt((n-arp)/(n-arp-1))

       oldb <- beta
       pp1 <- p + icent
       bsbeta <- matrix(rep(0,pp1^2),ncol=pp1)
       
       
       ind <- 1:n1
       ii <- sample(ind,1,replace=TRUE)
       allbeta <- c()
       rhostar <- c()
       MSEstar <- c()
       for(nbk in 1:nbs){
            ind2 <- sample(ind,n1,replace=TRUE)
            ystar <- y[ii:(ii+arp-1)]
            estar <- ehat[ind2]
            
            for(i in (arp+1):n){
                 ypart <- 0
                 for(k in 1:arp){ypart <- ypart + phi[k]*ystar[i-k]}
                 xpart <- 0
                 for(j in 1:pp1){
                      for(k in 1:arp){
                           if(k == 1){
                                xpart <- xpart + beta[j]*(xcopy[i,j]-phi[k]*xcopy[i-k,j])
                           } else {
                                xpart <- xpart - beta[j]*phi[k]*xcopy[i-k,j]
                           }
                      }
                 }
                 ystar[i] <- ypart + xpart + estar[i-arp]
            }
            
            avey <- mean(ystar)
            sigestar <- sd(estar*(n/(n-1)))
            MSEstar <- c(MSEstar, sigestar ^ 2)
            ystar <- ystar - avey
                    
            d2fit <- durbin2fit(ystar,x,phi,method=method,scores=scores)
            dum3 <- d2fit$beta
            d1 <- avey - t(xbar)%*%dum3
            dum3 <- c(d1,dum3)
            
            uhat <- y - xcopy %*% dum3
            dum4 = lagmat(uhat, arp)
            u.y <- dum4[,1]
            u.x <- dum4[,-1]
            ufit <- lm(u.y ~ u.x)
            rhotmp <- ufit$coef[-1]
            rhostar <- rbind(rhostar, rhotmp)
            
            allbeta <- rbind(allbeta, dum3)
            
            for(i in 1:pp1){
                 for(j in 1:pp1){
                     bsbeta[i,j] <- bsbeta[i,j]+(dum3[i]-oldb[i])*(dum3[j]-oldb[j])/sigestar^2
                 }
            }
       }

       bsbeta <- bsbeta / nbs
       
       list(betacov = bsbeta, allbeta = allbeta, rhostar = rhostar, MSEstar = MSEstar)
}
