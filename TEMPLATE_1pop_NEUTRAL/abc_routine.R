#abcmoments = function(target,sigma,mu,nsim,tol,abctable,rej=F,tothap,num1,num2,len,projection,maxss)  #original
abcmoments = function(target,sigma,mu,nsim,tol,abctable,rej=F,proj.mat=NULL,iswt) 

#mu is the mean vector from the cavity distribution - 
#sigma is the covariance matrix for the cavity distribution - 

#NB this code is inefficient because we can actually re-use an initially simulated ABC table, without resimulating
#for each locus, but this just gets things going. 
{

############################################################################################################################################




	#Most of what is below is quite generic - you do not need to use this code. What it is doing is the 
	#standard ABC calc - first it checks validy of summary stats (actually, above), then rescales them to have
	#the same variance (most ABC packages do this), then calculates distance from target, defines a 
	#cutoff, and then fits a regression model to try to get more precise estimation of the mean 
	#and covariance of the parameter values.
	#Alternatively you could just feed in sstable and params above to any abc package, and 
	#then use mean() and var() on the resulting posterior distribution matrix of parameter values 
	#rescale	
	#names(target) = names(sstable)
	
	nparams = length(mu)
	params = abctable$params
	sstable = abctable$sstable
	nmiss = abctable$nmiss

	# Projections	
	if(!is.null(proj.mat)){
		sstable = make_fp_proj(proj.mat,sstable)
		target = make_fp_proj(proj.mat,target)
		nss = length(target)
	}
	
	if(nss != nparams) stop("number of projected summary stats not equal to number of params")

	sdvec = apply(rbind(sstable,target),2,sd)
	mvec = apply(rbind(sstable,target),2,mean)
	sstable2 = sstable
	sstable2 = sweep(sstable2,2,mvec) #subtract column mean
	sstable2 = sweep(sstable2,2,sdvec,"/") #divide by column sd
	target2 = (target-mvec)/sdvec
		
# calc euclidean distance
	dst = sqrt(apply(sweep(sstable2,2,target2)^2,1,sum)) 
			#sweep() subtracts scaled target summ stats from the scaled abc table of summ stats
			#then apply() sums the squared values from sweep across the summary stats
			#then take the square root. 
 	tol = tol*nsim/(nsim - nmiss) #being pedantic - so the tolerance takes into account the weighted out vals
 	                                # i.e. a 1% tolerance means that we actually take the 1% nearest points 
 	                                # assuming that the weighted out vals have infinite distance. 
	ntol = ceiling(tol * nsim)
	if(ntol <= 0 )ntol = 1
	if(ntol > nsim)ntol = nsim
	thresh = sort(dst)[ntol]
	wt = dst <= thresh
	sampsize = sum(wt)
	regwt = iswt[wt]
    	regwt = (1-dst[wt]^2/thresh^2)*regwt
	wss = sum(regwt)^2/sum(regwt^2) #the effective sample size in the tolerance region

	
#################################################################################################
	#(10.01.2019) Now working on modifying the script to estimate marginal likelihood for each site. I will use kernel density 
	#estimation with the Epanechnikov kernel we are using for parameter estimation. 
	#So I think the kernel density estimate is 
	#	(sum of weights)/(number of points)/(kernel density volume)
	#In the case of e.g. a Gaussian kernel, the kernel is normalised, so this reduces to
	#	(sum of weights)/(number of points)
	#A further complication is that the points themselves have an importance weight. So I think this gives
	#	(sum of product of weights)/(sum of IS weights)/(kernel density volume)
	#I have checked this with toy simulation in R for both uniform and Gaussian kernel, and looks OK. 
	#
	# a further (minor, hopefully) complication is that the treatment (and R tests) above are for the univariate case,
	# may need to think more about multivariate.

	e.dens.wt = epanvol(dst[wt],thresh,nss) # This gives normalised k.d. weight for each point (remember to put in thresh^nss term)

	# from the treatment on multivariate density estimation we then have (assuming normalised kernel) that
	# density = (sum of product of weights)/(sum of IS weights)
	
	mv.dens = sum(e.dens.wt*iswt[wt])/sum(iswt)

	log.mv.dens = log(mv.dens)
	
	if(!is.finite(wss))wss = 0
	
	if(wss < 2*max(c(nss,nparams))){
		print(paste("ess problem for",target,": returning prior"))
		est.mu = mu
		est.cov = sigma
		est.precmat = chol2inv(chol(sigma))
		est.sigmainv.mu = as.numeric(est.mu%*%est.precmat)
 		#outlist = list(sigmainv.mu=est.sigmainv.mu,sigmainv = est.precmat,mu = est.mu,
 						#sigma = est.cov,ess = ess,pss=pss,params=params,wvec=wvec)	
 		#return(outlist)		
 		outlist = list(sigmainv.mu=est.sigmainv.mu,sigmainv = est.precmat,mu = est.mu,sigma = est.cov,abctable = abctable,wss=wss,tilt.cons = 0)	
                                                                                                                           #a bit uncertain about tilt.cons = 0 - revisit?	
 		return(outlist)		
	}
	
#################################################################################################
	
    #cumbersome code below is to make sure that we can do the prediction OK (probably easier just to use coeffs...)
    ss.df = data.frame(sstable2[wt,])
    pwt = as.matrix(params[wt,])
	xvar.names <- paste("v",as.character(c(1:nss)),sep="")
	names(ss.df) <- xvar.names

	fmla <- as.formula(paste("pwt ~ ", paste("as.matrix(",xvar.names,")", collapse= "+")))
	fmla2 <- as.formula(paste("pwt ~ 1")) #this is just for rejection; nicer to do it like this so we keep all the
	                                      #code for regression and rejection the same other than logical expression below
	                                      #so we can be confident we are comparing like with like
	if(rej){
		fit1 <- lm(fmla2,data=ss.df,weight=regwt)
		f.rank = 1 #number of columns in the design matrix; just fitting the intercept
		#print(paste("fmla2",ilocus))
	}else{
		fit1 <- lm(fmla,data=ss.df,weight=regwt)
		f.rank = nss + 1 #intercept plus covariates
		#print(paste("fmla",ilocus))
	}#	nparams = length(mu.cav)
#	nss = dim(target)[2]
	
	#if(nss != nparams) stop("number of projected summary stats not equal to number of params")
	
#	
#	if(wss <= nparams + 2 || wss <= f.rank){
#		print(paste("wss problem for",target,": returning prior"))
#		est.mu = mu
#		est.cov = sigma
#		est.precmat = chol2inv(chol(sigma))
#		est.sigmainv.mu = as.numeric(est.mu%*%est.precmat)
# 		outlist = list(sigmainv.mu=est.sigmainv.mu,sigmainv = est.precmat,mu = est.mu, sigma = est.cov,abctable = abctable)	
# 		return(outlist)			
#	}
#	
	
	#from is version, to check below: p.cent.left = sweep(p.cent,1,wvec,"*")
	#from is version est.cov = (t(p.cent.left) %*% p.cent)/sum(wvec)
	#calc above should be equivalent to est.cov = (t(p.cent) %*% diag(wvec) %*% p.cent)/sum(wvec)
	#need to do it this way for efficiency...
	
	leftcalc = sweep(fit1$residuals,1,regwt,"*")
	est.cov = (t(leftcalc) %*% fit1$residuals)/sum(regwt)
	
  	# est.cov = (t(fit1$residuals) %*% diag(regwt) %*% fit1$residuals)/sum(regwt)
    	est.cov = est.cov*wss/(wss - f.rank) 
    #I'm making the assumption that we can use wss in place of sampsize for the weighted case. This may need some
    #thought, and revisting, but it works OK for the importance sampling example I have worked on. 
  	#This is the estimated covariance matrix taking into account the regression of the parameter values on the summary stats.
  	#I.e. if we did not do regression then nss would be zero and it would the standard estimator but with wss in place of 
  	#sample size
 	#Ideally we want an unbiased estimate of the precision matrix, and this can be done as follows:

	est.precmat = chol2inv(chol(est.cov)) #find the inverse
	est.precmat = (wss-nparams - 2)/(wss-1)*est.precmat #this gives an unbiased estimate of precision matrix - again using wss in 
	                                                    #place of actual sample size...
	
	target.df = data.frame(matrix(target2,nrow=1,byrow=T))
	names(target.df) = xvar.names
	est.mu = predict(fit1,target.df)
 	       
 	est.sigmainv.mu = as.numeric(est.mu%*%est.precmat) 
 	#abctable$wss = wss
 	list(sigmainv.mu=est.sigmainv.mu,sigmainv = est.precmat,mu = est.mu,sigma = est.cov,abctable = abctable,wss=wss,tilt.cons = log.mv.dens,ntol=ntol) #Rita		


}









