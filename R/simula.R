simula <-
function(x,y,arp,nbs,nbscov, conf, method,scores=scores) {
    upper <- .99
    lower <- -.99
    n <- length(y)
    p <- length(x[1,])
    df <- n - p - arp
    icent <- 0
    if (p > 1) {
      icent <- 1
    }
    
    xcopy <- x
    
    ones <- rep(1,n)
    proj1 <- ones %*% t(ones) / n
    x2 <- x[,2:p]
    xbar <- apply(x2,2,mean)
    xc <- x[,2:p] - proj1 %*% x[,2:p]
    x <- xc
    p <- p - 1
    
    ##   durb1aa
    pcent <- p + icent
    
    dfit <- durbin1fit(y,xcopy,arp,method=method,scores=scores)
    adjphi <- dfit$coef[2:(arp + 1)]
    for (j in 1:arp) {
      if (adjphi[j] < lower) {
        adjphi[j] <- lower
      }
      if (adjphi[j] > upper) {
        adjphi[j] <- upper
      }
    }
    rho1 <- adjphi
    #     return(adjphi)
    
    ybar <- mean(y)
    yc <- y - ybar
    d2fit <- durbin2fit(yc,xc,adjphi,method=method,scores=scores)
    
    beta <- d2fit$beta
    b0 <- ybar - t(xbar) %*% beta
    
    allb <- c(b0,beta)
    
    
    ####
    cnt <- 0
    ic <- 0
    adjar <- adjphi
    #   beginning of bs loop
    while (ic == 0) {
      holdr <- adjar
      
      cnt <- cnt + 1
      np <- p + icent
      bsbias <- boot1(y,adjar,arp,nbs,xcopy,allb,method=method,scores=scores)
      
      for (k in 1:arp) {
        if (k == 1) {
          hild <- 0
        }
        adjar[k] <- adjphi[k] + bsbias[k]
        if (adjar[k] < lower) {
          adjar[k] <- lower
        }
        if (adjar[k] > upper) {
          adjar[k] <- upper
        }
      }
      diff <- adjar - holdr
      
      d2fit2 <- durbin2fit(yc,xc,adjar,method=method,scores=scores)
      beta2 <- d2fit2$beta
      b02 <- ybar - t(xbar) %*% beta2
      
      allb <- c(b02,beta2)
      check1 <- 0
      metric <- 0
      
      for (k in 1:arp) {
        if (diff[k] < 0) {
          check1 <- check1 + 1
        }
        metric <- metric + diff[k] ^ 2
      }
      metric <- sqrt(metric)
      if (check1 == arp) {
        adjar <- holdr
      }
      
      if (((metric <= .01) &&
           (abs(diff[1]) <= .01)) | (cnt > 8)) {
        ic <- 1
      }
    }
    
    ###########
    
    d2fit <- durbin2fit(yc,xc,adjar,method=method,scores=scores)
    
    beta <- d2fit$beta
    b0 <- ybar - t(xbar) %*% beta
    sigd2 <- d2fit$sigma
    mse <- c(sigd2 ^ 2)
    
    allb <- c(b0,beta)
    
    

    ##  bs cov mat
    bscov <- boot2(y,xcopy,adjar,allb,nbscov,method=method,scores=scores)
    betacov <- bscov$betacov * mse
    sesbeta <- diag(betacov) ^ (1 / 2)
    tees <- allb / sesbeta
    pvals <- 2 * (1 - pt(abs(tees),df))
    tabbeta <- cbind(allb,sesbeta,tees,pvals)
    colnames(tabbeta) <- c("beta","SE","t-ratio","p-value")
    rname <- c()
    for (j in 1:p) {
      rname <- c(rname,paste("beta_",j))
    }
    rownames(tabbeta) <- c("Intercept",rname)
    
    ### .99 flag ###
    flag99 <- 0
    if (any(adjar >= 0.99)) {
      flag99 <- 1
    }
    
    
    rhostar <- bscov$rhostar
    rhobias <- adjar - rho1
    k.multi <- (colMeans(rhostar) + rhobias) / colMeans(rhostar)
    # k.multi <- (rho1 + rhobias) / rho1
    rhostar2 = t(t(rhostar) + rhobias)
    rhostar <- t(t(rhostar) * k.multi) # multiply each row of rhostar matrix by the vector k.multi
    MSEstar <- bscov$MSEstar
    
    ### (a) ###
    rhocov1 <- matrix(rep(0,arp^2),ncol=arp)
    for(nbk in 1:nbscov){
      rhotmp <- rhostar[nbk,]
      for(i in 1:arp){
        for(j in 1:arp){
          rhocov1[i,j] <- rhocov1[i,j]+(rhotmp[i]-adjar[i])*(rhotmp[j]-adjar[j])
        }
      }
    }
    rhocov1 <- rhocov1 / nbscov
    serho1 <- diag(rhocov1) ^ (1 / 2)
    rho_CI_1 <- rbind(adjar - qt(1-(1-conf)/2, df) * serho1, adjar + qt(1-(1-conf)/2, df) * serho1)
    
    ### (b) ###
    rhocov2 <- matrix(rep(0,arp^2),ncol=arp)
    for(nbk in 1:nbscov){
      rhotmp <- rhostar[nbk,]
      MSEtmp <- MSEstar[nbk]
      for(i in 1:arp){
        for(j in 1:arp){
          rhocov2[i,j] <- rhocov2[i,j]+(rhotmp[i]-adjar[i])*(rhotmp[j]-adjar[j]) / MSEtmp
        }
      }
    }
    rhocov2 <- rhocov2 * mse / nbscov
    serho2 <- diag(rhocov2) ^ (1 / 2)
    rho_CI_2 <- rbind(adjar - qt((1-conf)/2, df) * serho2, adjar + qt(1-(1-conf)/2, df) * serho2)
    # ### (c) ###
    rho_CI_3 <- apply(rhostar2, 2, quantile, probs = c((1-conf)/2,1-(1-conf)/2))
    if (arp == 1){
      names(rho_CI_1) = c('LB','UB')
      names(rho_CI_2) = c('LB','UB')
      names(rho_CI_3) = c('LB','UB')
    } else {
      colnames(rho_CI_1) = c('LB','UB')
      colnames(rho_CI_2) = c('LB','UB')
      colnames(rho_CI_3) = c('LB','UB')
    }
    
    ### residuals and fitted.values
    ypart <- nurho(y, adjar)
    xpart <- wrho(xcopy, adjar)
    ehat <- ypart - xpart %*% allb
    fitted.values <- y[(arp+1):n] - ehat
   
   
    list(
      coefficients = unname(allb), rho1 = unname(rho1), adjar = unname(adjar), mse = mse, rho_CI_1 = rho_CI_1, rho_CI_2 = rho_CI_2, rho_CI_3 = rho_CI_3,
      betacov = betacov, tabbeta = tabbeta, flag99 = flag99, residuals = ehat, fitted.values = fitted.values
    )
  }
