hmmat <-
function(vecss,k){

      n <- sum(vecss)
      xmat <- matrix(rep(0,n*2*k),ncol=(2*k))
  
      ic <- 1
      ir <- 1

      for(i in 1:k){
           ni <- vecss[i]
           xmat[ir:(ir+ni-1),ic] <- 1
           if(i==1){
                 xmat[ir:(ir+ni-1),ic+1] <- 1:ni
           } else {
                 xmat[ir:(ir+ni-1),ic+1] <- 0:(ni-1)
           }
           if(i > 1){
                  ics <- 1
                  for(j in 1:(i-1)){
                        last <- xmat[ir-1,ics+1]
                        xmat[ir:(ir+ni-1),ics] <- 1
                        xmat[ir:(ir+ni-1),ics+1] <- last + (1:ni)
                        ics <- ics+2
                   }
            }
            ic <- ic + 2
            ir <- ir + ni
        }
        return(xmat)
}
