get_fp_proj = function(tolx=NULL,sstable,params,target=NULL)
{

	#June 2016. At this point we'll project the data and create a transformed sstable and target, but keep the same 
	#names so we don't have to modify the code further below. This has the disadvantage of using the data twice,but
	#keeps the coding easy, otherwise have to split data and then rename things, etc
	
	if(!is.null(tolx))target = as.numeric(target)
	nss = dim(sstable)[2]
	niter = dim(sstable)[1]
	if(!is.null(tolx)){
		scaling = numeric(nss)
		for(j in 1:nss){
			v1 = sstable[,j]
			v1 = abs(v1-target[j])
			v1 = sort(v1)
			scaling[j] = v1[niter*tolx]
		}
		scaling = ifelse(scaling == 0,1,scaling)
		centering = target
		sstable2 = sstable
		sstable2 = sweep(sstable2,2,centering) #subtract target value
		sstable2 = sweep(sstable2,2,scaling,"/") #divide by scaling
		target2 = (target-centering)/scaling #so this should be 0
	}else{
		scaling = apply(sstable,2,mad)
		scaling = ifelse(scaling == 0,1,scaling)
		centering = apply(sstable,2,median)
		sstable2 = sstable
		sstable2 = sweep(sstable2,2,centering) #subtract column median
		sstable2 = sweep(sstable2,2,scaling,"/") #divide by column mad	
	}
		
 	#calc euclidean distance
	if(!is.null(tolx))dst = sqrt(apply(sweep(sstable2,2,target2)^2,1,sum)) 
			#sweep() subtracts scaled target summ stats from the scaled abc table of summ stats
			#then apply() sums the squared values from sweep across the summary stats
			#then take the square root. 
	if(is.null(tolx)){
		wt = rep(T,niter)
		sampsize = niter
	}else{
		ntol = ceiling(tolx * niter)
		if(ntol <= 0 )ntol = 1
		if(ntol > niter)ntol = niter
		thresh = sort(dst)[ntol]
		wt = dst <= thresh
		sampsize = sum(wt)
	}
	
    #cumbersome code below is to make sure that we can do the prediction OK (probably easier just to use coeffs...)
    ss.df = data.frame(sstable2[wt,])
    pwt = params[wt,]
	xvar.names <- paste("v",as.character(c(1:nss)),sep="")
	names(ss.df) <- xvar.names

	fmla <- as.formula(paste("pwt ~ ", paste(xvar.names, collapse= "+")))
	fit1 <- lm(fmla,data=ss.df)
	#fitted(fit1) should now contain the projected summary stats
	coeff.ok = ifelse(is.finite(fit1$coeff),fit1$coeff,0)# we'll just zero out problematic columns
	list(coeff = coeff.ok,centering = centering, scaling = scaling)
}

make_fp_proj = function(trans,target)
{
	scaling = trans$scaling
	centering = trans$centering
	m1 = trans$coeff
	if(is.null(dim(m1)))stop("projection matrix not matrix")
	if(is.null(dim(target))){
		if(length(target) != dim(m1)[1] - 1)stop("vector argument has incompatible length")
		target = (target-centering)/scaling #so this should be 0		
		return(as.vector(t(c(1,target))%*%m1))
	}
	if(dim(target)[2] + 1 != dim(m1)[1])stop("matrix argument has incompatible number of columns")
	sstable2 = target
	sstable2 = sweep(sstable2,2,centering) #subtract target value
	sstable2 = sweep(sstable2,2,scaling,"/") #divide by scaling
	return(cbind(1,sstable2)%*%m1)
	
}
	
