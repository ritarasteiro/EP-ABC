##############################################################################################################
##NEUTRAL SCENARIOS
#This script uses the posteriores calculated in epABC and uses them as priors. 
##It creates 2 parameter files: (i) from cavity distribution to be used in the ep-abc, (ii) with times converted to be used in SLIM
##############################################################################################################
args = commandArgs(trailingOnly=TRUE)
options(scipen=999) # convert scientific notation to numeric (to work both in slim and msprime)
library(mnormt)

Npara=13
Nloci = as.numeric(args[4])

# Number of simulations
numsim = as.numeric(args[1]) #number of ABC simulations
SCALE=as.numeric(args[2])
#---------------------
# Fixed priors
#---------------------
#generation scale parameter
GEN=as.numeric(args[3])

#time in generations
# From Boitard: Modifying T would allow population size changes to occur on a longer or shorter period in the past, 
#and modifying a would allow to describe more precisely one specific part of the history, playing on the ratio
# between the length of recent versus old time windows.

 T=50000/GEN
 a=0.1
 L=12 #nb windows
 t=array(0)
 for (i in 1:L-1){
 t[i]=ceiling((exp(log(1 + a*T)*i/(L - 1)) - 1)/a)
 }
#[1]    10    29    66   139   282   563  1114  2195  4318  8486 16667


#------------------------------------------
# Priors drawn from the cavity distribution
#---------------------
#we build up the posteriors by initially starting with the prior then 
#getting the posterior for one locus, then using this as the prior for the next locus
#etc. We then refine the posterior in the subsequent sweeps through the data. 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("covfunc.R")
#time points (centroids)
time=array(0)
time[1]=t[1]/2
for(i in 2:11){
time[i]=(t[i]-t[i-1])/2
}
time[12]=((100000/GEN))

pos=log(time) #time on log
#Use this covariance matrix to give a prior to N
cov1 = covfunc(pos,1.5,2.0,0.1) #CHANGED lower lambda=1.5
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#set priors. 
#Prior: N0-N10,Nanc,REC
#sd
Prior.matrix = diag(c(rep(0,12),0.5)) 
Prior.matrix[1:12,1:12]=cov1 #impose the covariance prior for N
Prior.sigmainv=chol2inv(chol(Prior.matrix)) #precision matrix (Prior.matrix inverted)

#means
Priors=log(c(rep(30000,12),2.3e-8)) 
Prior.sigmainv.mu = as.vector(Priors%*%Prior.sigmainv) #centre it around nothing very interesting

#first we set the natural parameters to values that give a guess at the posterior. 
# I've set these guesses below and then we need to work backwards to get the natural parameters
guess.posterior.cov = diag(rep(0.2,Npara))# increase from 0.1 to 0.2 

#guess.posterior.mean = log(c(rep(20000,12),4.5e-8))
############################
target=as.matrix(read.table(paste(args[5],"_targetSS.txt",sep="")))


	params.proj <-matrix(scan("../PROJECTION/params.proj.txt"),ncol=Npara,byrow=T) 
	sstable.proj = as.matrix(read.table("../PROJECTION/SS_proj.txt")) #summary stats
	numsim.proj=dim(params.proj)[1]
	
	gwt.proj = apply(is.finite(sstable.proj),1,all) 
	nmiss.proj = numsim.proj - sum(gwt.proj)
	if(nmiss.proj > 0){
		sstable.proj = sstable.proj[gwt.proj,]
		params.proj = params.proj[gwt.proj,]
	}


nmiss=nmiss.proj	
tol=0.05
nsim=numsim.proj
sstable=sstable.proj
params=params.proj
est.mu.all=matrix(nrow=Nloci,ncol=Npara)

for (i in 1:Nloci){

dst = sqrt(apply(sweep(sstable.proj,2,target[i,])^2,1,sum)) 
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
#	regwt = iswt[wt]
 #   	regwt = (1-dst[wt]^2/thresh^2)*regwt
#	wss = sum(regwt)^2/sum(regwt^2)
	
	
	 #cumbersome code below is to make sure that we can do the prediction OK (probably easier just to use coeffs...)
    ss.df = data.frame(sstable)
    pwt = as.matrix(params)
	nss=length(target[i,])
	xvar.names <- paste("v",as.character(c(1:nss)),sep="")
	names(ss.df) <- xvar.names

	fmla <- as.formula(paste("pwt ~ ", paste("as.matrix(",xvar.names,")", collapse= "+")))
		fit1 <- lm(fmla,data=ss.df) #,weight=regwt)
		f.rank = nss + 1 #intercept plus covariates
		#print(paste("fmla",ilocus))
	
	target.df = data.frame(matrix(target[i,],nrow=1,byrow=T))
	names(target.df) = xvar.names
	est.mu = predict(fit1,target.df)
	est.mu.all[i,]=est.mu
 }
est.mu.all.mean=colMeans(est.mu.all)

guess.posterior.mean = est.mu.all.mean

################################################################################
#guess.posterior.mean = as.vector(t(read.table("mean_start_values.txt", sep="")))
#####
SigmaInv.sum = chol2inv(chol(guess.posterior.cov)) - Prior.sigmainv 
#this is going to be the inverse of the covariance matrix of the posterior
#(minus the prior) for the parameters


SigmaInv.mu.sum = chol2inv(chol(guess.posterior.cov))%*%guess.posterior.mean - Prior.sigmainv.mu 
			       #this is going to be the mean for posterior distribution (minus the prior) of 
                               #the parameters divided by the covariance matrix.

SigmaInv.mu.cav = (Nloci-1)/Nloci*SigmaInv.mu.sum + Prior.sigmainv.mu
SigmaInv.cav = (Nloci-1)/Nloci*SigmaInv.sum + Prior.sigmainv
sigma.cav = chol2inv(chol(SigmaInv.cav))
mu.cav = sigma.cav%*%SigmaInv.mu.cav 


#simulate from cavity distribution (acts as a prior)
params_abctable = rmnorm(numsim, mean = as.numeric(mu.cav), varcov = sigma.cav) 

#gives the denominator in the IS calculation
dens.params_abctable = dmnorm(params_abctable,mean=as.numeric(mu.cav), varcov = sigma.cav) 

#bring back from natural logs
ep_prior = exp(params_abctable)

#---------------------
# Priors for simulations
#---------------------
# Population sizes
N=ceiling(ep_prior[,1:11]/SCALE)

#Ancestral population
Nanc=ceiling(ep_prior[,12]/SCALE)


#################################################################
#MUTATION, RECOMBINATION AND SELECTION
#################################################################
#recombination parameter (relative to mutation rate) 
mut=3.5e-9*SCALE  # nucleotide per generation (cichlids)
rec <-1/2*(1-(1-(2*ep_prior[,13]))^SCALE)# using formula suggested in SLIM manual
#---------------------
# Matrix of parameters for each simulation
#---------------------
params.sims <- data.frame(
	N0 = N[,1],
	N1 = N[,2],
	N2 = N[,3],
	N3 = N[,4],
	N4 = N[,5],
	N5 = N[,6],
	N6 = N[,7],
	N7 = N[,8],
	N8 = N[,9],
	N9 = N[,10],
	N10 = N[,11],
	Nanc = Nanc,
	T1=ceiling(t[1]/SCALE),
	T2=ceiling(t[2]/SCALE),
	T3=ceiling(t[3]/SCALE),
	T4=ceiling(t[4]/SCALE),
	T5=ceiling(t[5]/SCALE),
	T6=ceiling(t[6]/SCALE),
	T7=ceiling(t[7]/SCALE),
	T8=ceiling(t[8]/SCALE),
	T9=ceiling(t[9]/SCALE),
	T10=ceiling(t[10]/SCALE),
	T11=ceiling(t[11]/SCALE),
	MUT=mut,
	REC=rec)# using formula suggested in SLIM manual


dir.create("iter1")

#For ep_abc
write.table (params_abctable, "iter1/params_abctable.txt", col.names=F,row.names=F)
write.table (dens.params_abctable, "iter1/dens.params_abctable.txt", col.names=F,row.names=F)
write.table (sigma.cav, "iter1/sigma.cav", col.names=F,row.names=F)
write.table (mu.cav, "iter1/mu.cav", col.names=F,row.names=F)
write.table (SigmaInv.sum,"iter1/SigmaInv.sum", col.names=F,row.names=F)
write.table (SigmaInv.mu.sum, "iter1/SigmaInv.mu.sum", col.names=F,row.names=F)
#save(c(params_abctable,dens.params_abctable,sigma.cav,mu.cav,SigmaInv.sum,SigmaInv.mu.sum), file="iter1/iter1_ep.RData")

# Prior
save(Prior.sigmainv,Prior.sigmainv.mu, file="iter1/Prior.RData")
#write.table (Prior.sigmainv,"iter1/Prior.sigmainv", col.names=F,row.names=F)
#write.table (Prior.sigmainv.mu, "iter1/Prior.sigmainv.mu", col.names=F,row.names=F)

#For sims
write.table (params.sims, "iter1/params.sims1.txt", col.names=F,row.names=F)
write(colnames(params.sims),"params.colnames.sims.txt", ncol=dim(params.sims)[2])

iter=1
write(iter,"iter_nb")
write(iter, paste("iter",iter,"/iter_nb",sep=""))

#write file with values sampled from prior
write.table(guess.posterior.mean,"mean_start_values.txt")




