nurho <-
function(yc,adjphi){

     arp <- length(adjphi)
     lagdata <- lagmat(yc,arp)
     yresp <- lagdata[,1]
     ylag <- lagdata[,2:(arp+1)]
     ylag <- as.matrix(ylag)
     nurho <- yresp - ylag%*%adjphi
     return(nurho)
}
