simulacorrection <-
function(x,y,arp,nbs,nbscov,method,scores) {
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
    
    ##   Durbin stage 1
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
    
    ## Durbin stage 2
    ybar <- mean(y)
    yc <- y - ybar
    d2fit <- durbin2fit(yc,xc,adjphi,method=method,scores=scores)
    
    beta <- d2fit$beta
    b0 <- ybar - t(xbar) %*% beta
    
    allb <- c(b0,beta)
    
    adjar <- adjphi
    rho1 <- adjphi
    holdr <- adjar
    
    ## only perform bootstrap once
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
    
    ### sse ###
    ypart <- nurho(y,adjar)
    xpart <- wrho(xcopy,adjar)
    ehat <- ypart - xpart %*% allb
    sse <- sum(ehat ^ 2)
    ehat2 <- y - xcopy %*% allb
    sse2 <- sum(ehat2 ^ 2)
    ### CI ###
    
    cv <- abs(qnorm((1 - 0.95) / 2))
    
    CI_adj_rho <- rhoci2(n,adjar,cv)
    adj_mid <- (CI_adj_rho[2] + CI_adj_rho[3]) / 2
    CI_init_rho <- rhoci2(n,rho1,cv)
    init_mid <- (CI_init_rho[2] + CI_init_rho[3]) / 2
    if (adj_mid <= 0.95) {
      adjar <- adj_mid
    } else {
      adjar <- init_mid
    }
    
   
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
    pvals <- 2 * (1 - pt(abs(tees), df))
    tabbeta <- cbind(allb, sesbeta, tees, pvals)
    colnames(tabbeta) <- c("beta","SE","t-ratio","p-value")
    rname <- c()
    for (j in 1:p) {
      rname <- c(rname, paste("beta_", j))
    }
    rownames(tabbeta) <- c("Intercept", rname)
    flag99 <- 1
    
    ypart <- nurho(y, adjar)
    xpart <- wrho(xcopy, adjar)
    ehat <- ypart - xpart %*% allb
    fitted.values <- y[2:n] - ehat
    rho_CI_1 <- c(NA,NA)
    rho_CI_2 <- c(NA,NA)
    rho_CI_3 <- c(NA,NA)
    list(
      coefficients = unname(allb), rho1 = unname(rho1), adjar = unname(adjar), mse = mse, rho_CI_1 = rho_CI_1, rho_CI_2 = rho_CI_2,
      rho_CI_3 = rho_CI_3, betacov = betacov, tabbeta = tabbeta, flag99 = flag99, residuals = ehat, fitted.values = fitted.values
    )
  }
