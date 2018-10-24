wrho <-
function(xc,adjphi){

     arp <- length(adjphi)
     n <- length(xc[,1])

     s1 <- arp+1
     s2 <- n

     x1 <- lagx(xc,s1,s2)

     for(j in 1:arp){
         s1 <- s1 - 1
         s2 <- s2 - 1
         xt <- lagx(xc,s1,s2)*adjphi[j]
         x1 <- x1 - xt
     }
     wrho <- x1
     return(wrho)
}
