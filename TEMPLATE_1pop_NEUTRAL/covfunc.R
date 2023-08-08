tridiag <-
function(upper, lower, main){
    out <- matrix(0,length(main),length(main))
    diag(out) <- main
    indx <- seq.int(length(upper))
    out[cbind(indx+1,indx)] <- lower
    out[cbind(indx,indx+1)] <- upper
    return(out)
}

covfunc <- 
#pos is just the time points (on a log 10 scale, probably best) 
#don't need to worry too much about units of pos because lambda scales it
#lambda is some scaling value - big numbers make changes smoother, small
#numbers make it wigglier
#On the diagonal the variance is sig1 + sig2
#Off the diagonal the covariance declines with distance between time points, 
#modulated by sig1 and lambda
#As lambda -> infinity, then you will have a diagonal matrix with variance sig2
function(pos,lambda,sig1,sig2){
	n1 = length(pos)
	cmat = matrix(nrow=n1,ncol=n1)
	for(j in 1:n1){
		for(i in 1:n1){
			cmat[i,j] = sig1^2*exp(-1/(2*lambda^2)*(pos[i] - pos[j])^2) + ifelse(i==j,1,0)*sig2^2
		}
	}
	return(cmat)
}
