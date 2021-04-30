durbin1xy <-
function(y,x,arp){

      lagy <- lagmat(y,arp)
      y2 <- lagy[,1]
      part1 <- lagy[,2:(arp+1)]
      n <- length(y)
      p <- length(x[1,])

      s1 <- arp + 1
      s2 <- n
      part2 <- lagx(as.matrix(x[,2:p]),s1,s2)
      allx <- cbind(part1,part2)

      for(j in 1:arp){
         s1 <- s1 -1
         s2 <- s2 -1
         part3 <- lagx(as.matrix(x[,2:p]),s1,s2)
         allx <- cbind(allx,part3)
      }

      durbin1xy <- cbind(allx,y2)
      return(durbin1xy)
}
