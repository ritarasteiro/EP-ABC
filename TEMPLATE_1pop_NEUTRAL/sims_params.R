##############################################################################################################
##This script uses the posteriores calculated in epABC and uses them as priors. 
##It creates 2 parameter files: (i) from cavity distribution to be used in the ep-abc, (ii)  with times converted to be used in SLIM and msprime
##############################################################################################################
args = commandArgs(trailingOnly=TRUE)
create.param.sims = function(iter,params_abctable,SigmaInv.sum,SigmaInv.mu.sum,numsim,Npara,SCALE,GEN) 

{
options(scipen=999) # convert scientific notation to numeric (to work both in slim and msprime)

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

ep_prior = exp(params_abctable)
#---------------------
# Priors for simulations
#---------------------
# Population sizes
N=ceiling(ep_prior[,1:11]/SCALE)

#Ancestral population
Nanc=ceiling(ep_prior[,12]/SCALE)

##################################################################
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


#For ep_abc
write.table (SigmaInv.sum,paste("iter",iter,"/SigmaInv.sum",iter,sep=""), col.names=F,row.names=F)
write.table (SigmaInv.mu.sum, paste("iter",iter,"/SigmaInv.mu.sum",iter,sep=""), col.names=F,row.names=F)
#
write(iter,"iter_nb")
write(iter, paste("iter",iter,"/iter_nb",sep=""))
#For sims
write.table (params.sims, paste("iter",iter,"/params.sims",iter,".txt",sep=""), col.names=F,row.names=F)

}

