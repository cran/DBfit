fullr <- function(x,p1){

#   This code assumes that the first p1 columns of x are lineraly independent
#   This should be checked before this function is called.

	p <- dim(x)[2]
	x2 <- x[,1:p1]
	rx <- qr(x)$rank

	if(rx==p){
		x2 <- x
	} else {
		j <- p1 + 1
		r2 <- p1
		while(j <= p){
			xtmp <- cbind(x2,x[,j])
			rtmp <- qr(xtmp)$rank
			if(rtmp > r2){
				r2 <- r2 + 1
				x2 <- xtmp
			}
			j <- j+1
		}
	}
	return(x2)
}
