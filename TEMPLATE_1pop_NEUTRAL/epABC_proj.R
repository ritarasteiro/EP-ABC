#CHANGE 20 Nov 2019
# Using projection functions in fp_proj.R
# It is parallel ABC, but it is still in a Loop: TODO parallelise

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#This script will run EP-ABC, following Barthelme and Chopin.
#Because the possible values for the parameters are bounded by zero, and since 
#everything is assuming a Gaussian distribution (posteriors and priors) it seems 
#natural to work on the log scale (natural logs rather than log10). 

#It uses the 'cavity distribution' that can be thought of as the estimate of the prior for a particular locus,
#when using all the other data. 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#LOAD libraries and source code
library(mnormt)
source("abc_routine.R")
source("marg_lik.R")
source("fp_proj.R")

##READ the arguments listed at the command line
args = commandArgs(trailingOnly=TRUE)

#TARGET
target=as.matrix(read.table(paste(args[1],"_targetSS.txt",sep="")))


#INITIALIZE THINGS:
Niter = as.numeric(args[2]) #this is the number of sweeps through the loci 
Nloci = dim(target)[1]
Smooth =Nloci*5
numsim = as.numeric(args[3]) #number of ABC simulations
tol = as.numeric(args[4]) #ABC we'll choose the ?% nearest points

#PRIORS

load("iter1/Prior.RData")
#Prior.sigmainv=read.table("iter1/Prior.sigmainv")
#Prior.sigmainv.mu=scan("iter1/Prior.sigmainv.mu")

#NUMBER of PARAMETERS
Npara = length(Prior.sigmainv.mu)

#CREATE ABC TABLE

#we build up the posteriors by initially starting with the prior then 
#getting the posteriors for all locus, join them together to create a global posterior and then using this as the prior for the iteration
#etc. We then refine the posterior in the subsequent sweeps through the data. 
abctable = NULL

#moved here: part of abc moments where it creates new abctable
# Now we are re-using an initially simulated ABC table, without resimulating for each locus.
iter_nb=scan("iter_nb")

#for iteration 1	
if (iter_nb == 1) {
	params <-matrix(scan("iter1/params_abctable.txt"),ncol=Npara,byrow=T) 
	dens.params = scan("iter1/dens.params_abctable.txt")#gives the denominator in the IS calculation
	sstable = as.matrix(read.table("iter1/SS1.txt")) #summary stats
	SigmaInv.sum=matrix(scan("iter1/SigmaInv.sum"),ncol=Npara,byrow=T)
	SigmaInv.mu.sum=matrix(scan("iter1/SigmaInv.mu.sum"),nrow=Npara,byrow=T)
	
	## this just checks that we have valid ss
	gwt = apply(is.finite(sstable),1,all) 
	nmiss = numsim - sum(gwt)
	if(nmiss > 0){
		sstable = sstable[gwt,]
		params = params[gwt,]
		dens.params = dens.params[gwt]
	}
	#PROJECT the data and create a transformed sstable and target: create a list that will be used by make_fp_poj in abc_routine.R
	# from here on all iterations will always be using the same projection
	params.proj <-matrix(scan("../PROJECTION/params.proj.txt"),ncol=Npara,byrow=T) 
	sstable.proj = as.matrix(read.table("../PROJECTION/SS_proj.txt")) #summary stats
	numsim.proj=dim(params.proj)[1]
	
	gwt.proj = apply(is.finite(sstable.proj),1,all) 
	nmiss.proj = numsim.proj - sum(gwt.proj)
	if(nmiss.proj > 0){
		sstable.proj = sstable.proj[gwt.proj,]
		params.proj = params.proj[gwt.proj,]
	}

	proj.mat = list()
	for(j in 1:Nloci){
		proj.mat[[j]] = get_fp_proj(tolx = 0.05,sstable=sstable.proj,params=params.proj,target=target[j,])
		}
	#save projection list
	save(proj.mat, file="proj.mat.RData")

#for iterations > 1
} else{
	params <-matrix(scan(paste("iter",iter_nb,"/abc_params",iter_nb,".txt",sep="")),ncol=Npara,byrow=T) 
	dens.params = scan(paste("iter",iter_nb,"/abc_dens.params",iter_nb,".txt",sep=""))#gives the denominator in the IS calculation
	sstable = as.matrix(read.table(paste("iter",iter_nb,"/SS",iter_nb,".txt",sep=""))) #summary stats
	SigmaInv.sum=matrix(scan(paste("iter",iter_nb,"/SigmaInv.sum",iter_nb,sep="")),ncol=Npara,byrow=T)
	SigmaInv.mu.sum=matrix(scan(paste("iter",iter_nb,"/SigmaInv.mu.sum",iter_nb,sep="")),nrow=Npara,byrow=T)
	
	## this just checks that we have valid ss
	gwt = apply(is.finite(sstable),1,all) 
	nmiss = numsim - sum(gwt)
	if(nmiss > 0){
		sstable = sstable[gwt,]
		params = params[gwt,]
		dens.params = dens.params[gwt]
	}
	#load projection list saved in iteration 1
	load("proj.mat.RData")
}

iswt = rep(1.0,numsim-nmiss)
ess = numsim - nmiss #effective sample size equal to sample size
abctable = list(params=params,dens.params = dens.params, sstable = sstable, nmiss = nmiss, ess = ess,wss = NA)


#START EP
for(iter in iter_nb:Niter) {
	print(paste("ITERATION",iter))
	SigmaInv.mu.cav = (Nloci-1)/Nloci*SigmaInv.mu.sum + Prior.sigmainv.mu
	SigmaInv.cav = (Nloci-1)/Nloci*SigmaInv.sum + Prior.sigmainv
	sigma.cav = chol2inv(chol(SigmaInv.cav))
	mu.cav = sigma.cav%*%SigmaInv.mu.cav 

	startagain = F # To make things quicker, as in Barthelme and Chopin, we re-use the ABC table, using importance sampling. 
	               #The basic idea is that, even though initially the parameters were simulated from some distribution with a particular 
	               #mean and variance, we can 'pretend' the parameter value were simulated under some other distribution as long
	               #as we know the importance weight. This is the probability density of the original points under the new model divided
	               #by the probability density of the original points under the original model they were simulated from. So, if the value is 
	               #0.5 for a particular simulated theta and rho you are saying that this would have been simulated half as often under
	               #the new model. So we can then use these importance weights to compute weighted means and covariances, but just based
	               #on an older set of data. However at some stage the distribution of weights gets very skewed, so that the effective
	               #sample size is very small, making the ABC unstable. In this case we recompute the table, and startagain is set to
	               #T
	
	if(!is.null(abctable)){ #if we have previously generated an abctable (i.e. we are in the middle of the simulation)
		iswt = dmnorm(params,mean=as.numeric(mu.cav),varcov=sigma.cav)/dens.params #this is the crucial bit that allows us to re-use the abc table
		ess = sum(iswt,na.rm=T)^2/sum(iswt^2,na.rm=T) #this computes an effective sample size; if all weights equal, it is the same as the size of the table
		print(paste("ess is",ess, ", iter nb",iter))		
		if(ess/(numsim-abctable$nmiss) < as.numeric(args[5])){ 
		startagain = T
		print(paste("startagain TRUE",ess)) #this might get annoying; I just have it to monitor things
		} else {
			startagain = F
			abctable$ess = ess
			print(paste("startagain FALSE",ess)) #this might get annoying; I just have it to monitor things
		}
	
	}
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#WHEN IT IS NEEDED TO SIMULATE AGAIN: startagain = T

	if(is.null(abctable) || startagain){#when this function is first used, abctable is null, so we start here first
	                                    #also if we start again
		
		sigma.trans=sigma.cav# to deal with extreme covariances
		#diag(sigma.trans)=1.4*diag(sigma.trans)
		sigma.trans=1.4*sigma.trans # use this one instead
		
		#simulate from cavity distribution (acts as a prior)
		params = rmnorm(numsim, mean = as.numeric(mu.cav), varcov = sigma.trans) 		
		dens.params = dmnorm(params,mean=as.numeric(mu.cav),varcov = sigma.trans) #gives the denominator in the IS calculation
		
		# create files and folders needed to run the simulators		
		
		if (iter==iter_nb){
		# create a folder
		dir.create(paste("iter",iter+1,sep=""))
		write.table(params,file=paste("iter",iter+1,"/abc_params",iter+1,".txt",sep=""), col.names=F,row.names=F)
		write.table(dens.params,file=paste("iter",iter+1,"/abc_dens.params",iter+1,".txt",sep=""), col.names=F,row.names=F)

		source("sims_params.R")
		create.param.sims(iter+1,params,SigmaInv.sum,SigmaInv.mu.sum,numsim,Npara,as.numeric(args[6]),as.numeric(args[7]))  
	
		} else{
		# create a folder
		dir.create(paste("iter",iter,sep=""))
		write.table(params,file=paste("iter",iter,"/abc_params",iter,".txt",sep=""),col.names=F,row.names=F)
		write.table(dens.params,file=paste("iter",iter,"/abc_dens.params",iter,".txt",sep=""),col.names=F,row.names=F)

		source("sims_params.R")
		create.param.sims(iter,params,SigmaInv.sum,SigmaInv.mu.sum,numsim,Npara,as.numeric(args[6]),as.numeric(args[7]))  
		}		
	print(paste("New input created"))
	stop("got to here")
	}
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	#A SWEEP THROUGH THE LOCI
	lambda1.sum = 0
	lambda2.sum = 0
	ess.av=0
	Norm.Cons = 0

	for(ilocus in 1:Nloci){ #goes locus by locus
		# Compute the parameters of the cavity distribution. This
		# is not standard for EP - I am using what Barthelme and Dehaene term average EP. 
		# In normal EP, each locus has its own posterior approximation - the overall posterior distribution is 
		# given as the product of a set of multivariate normals, one for each locus. The parameters of this
		# distribution is simply the sum of the natural parameters for each locus-specific distribution. So
		# in standard EP we create the cavity distribution for a particular locus by taking the natural parameters
		# for the current approximation to the posterior (without the prior) (i.e. SigmaInv.mu.sum and SigmaInv.sum) and subtract out
		# the natural parameter values for the current locus. Below, instead, what we are implicitly assuming is
		# that each locus contributes equally to the posterior, so we can obtain the natural parameter for the 
		# cavity distribution by multiplying the current posterior (without prior) by (Nloci -1)/Nloci. This has significant advantages
		# when we are getting noisy estimates for the mean and covariance because otherwise sometimes we can get negative
		# natural parameters (not allowed) for the locus specific posterior approximation. 
		# BTW The '(without the prior)' term in the explanations refers to the fact that I keep the prior natural paramters
		# separate in the average-EP calculation, and just put them in when we are computing the cavity distribution and 
		# in presenting the final results. 
		#..SigmaInv.mu.cav = (Nloci-1)/Nloci*SigmaInv.mu.sum + Prior.sigmainv.mu
		#..SigmaInv.cav = (Nloci-1)/Nloci*SigmaInv.sum + Prior.sigmainv
		
		#use simulation to get moments of the tilted distribution need to do this because the multivariate normal simulator
		# needs the standard parameters rather than the natural parameters. (Of course, could do this in the abc function)
		#..sigma.cav = chol2inv(chol(SigmaInv.cav))
		
		
		#..mu.cav = sigma.cav%*%SigmaInv.mu.cav 

		print(paste("LOCUS",ilocus))
		simres = abcmoments(target[ilocus,],sigma.cav,as.numeric(mu.cav),numsim,tol,abctable,rej=F,proj.mat=proj.mat[[ilocus]],iswt) 
		                                                 #this is because rmnorm() 
		                                                 #doesn't like a matrix mu.cav
		#abctable = simres$abctable #so the key output of abcmoments is a refined posterior distribution, returned just simply
		                           #as revised values for the natural parameters of the posterior # NOT USED: always using the same abctable (parallel)in each iteration
                       
		                         
		lambda1.sum = lambda1.sum + simres$sigmainv.mu - SigmaInv.mu.cav #this gives the approximating likelihood for each locus
		lambda2.sum = lambda2.sum + simres$sigmainv - SigmaInv.cav
		
		print(paste("wss",simres$wss))
		ess.av = ess.av + simres$wss #RITA
		#see Evernote discussion. Finally ended up with this. The Seeger review is most reliable (not Gelman review).
		Norm.Cons = Norm.Cons + simres$tilt.cons - phicalc(SigmaInv.mu.sum+Prior.sigmainv.mu,as.matrix(SigmaInv.sum+Prior.sigmainv)) + 
								phicalc(SigmaInv.mu.cav,as.matrix(SigmaInv.cav)) 
		
		#write(cbind(iter,ilocus,simres$ntol,simres$wss),file="ABCreport.txt",ncol=4,append=T)		
	}
	report_Norm.Cons = Norm.Cons + phicalc(SigmaInv.mu.sum+Prior.sigmainv.mu,as.matrix(SigmaInv.sum+Prior.sigmainv)) - phicalc(Prior.sigmainv.mu,as.matrix(Prior.sigmainv))

	ess.av = ess.av/Nloci
	print(paste("average ess is",ess.av, ", iter nb",iter))

	#GLOBAL POSTERIORS
	SigmaInv.mu.sum = (Smooth - 1)/Smooth*SigmaInv.mu.sum + lambda1.sum/Smooth
	SigmaInv.sum = (Smooth - 1)/Smooth*SigmaInv.sum + lambda2.sum/Smooth

#	REPORT RESULTS FOR THIS SWEEP 
	sigmat = chol2inv(chol(SigmaInv.sum + Prior.sigmainv)) #note now put prior back in for reporting the results of this sweep
	muvec = as.numeric(sigmat%*%(SigmaInv.mu.sum + Prior.sigmainv.mu)) #convert back to standard parameters (mean and covariance)
	sdev = sqrt(diag(sigmat)) #standard deviation
	smat = diag(1/sdev)
	corrmat = smat%*%sigmat%*%smat #gets correlation matrix rather than covariance (easier to visualise)
	index1 = upper.tri(corrmat)
	
	# output contains the marginal likelihood, THE vector of posterior means, the vector
	# of posterior standard deviations, and the correlation coefficients
	output = c(report_Norm.Cons,muvec,sdev,corrmat[index1]) #just reports this as a row of output (note the correlation matrix is flattened)
	write(output,file=paste("output",iter,".txt",sep=""),ncol=length(output),append=F)
}

if (iter==Niter){
write(iter,"iter_nb")
rm(list=ls())

}else{
rm(list=ls())
}
	

