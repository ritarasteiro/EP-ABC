# Changed to include new sumstats 12 Nov 2019

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#to give an argument in commad line
# R functions associated with calculating summary statistics
source("../files_to_run/Rfuncs.r")
library(Rcpp)
##library(Rfast)
Rcpp::sourceCpp("../files_to_run/diseqSS.cpp",cacheDir=paste(args[4],"/",sep="")) #to avoid compiling at every sim
#Rcpp::sourceCpp("../files_to_run/diseqSS.cpp")
library(mnormt)

samps<-as.integer(args[2])    
nsamp= length(samps)
nind=samps/2

# Read the MS output from SLiM
txt <- scan(args[1], what = character(0),sep = "\n", quiet = TRUE)

NBASES=as.integer(args[3]) #TODO

test <- read.ms.output_compiled(txt, nsam=samps)
#---------------------
# Calculate summary statistics
#---------------------
# Create matrix combining subpopulations

geno.all<- matrix(,ncol=length(unique(unlist(test$positions))[!is.na(unique(unlist(test$positions)))]),)

i=1
geno<-matrix(,length(unique(unlist(test$positions))[!is.na(unique(unlist(test$positions)))]), nrow=test$nsam[i])
colnames(geno)<-sort(unique(unlist(test$positions)))

if(length(test$gametes[[i]]) > 0)
geno[,as.character(test$positions[[i]])]<-test$gametes[[i]]
assign(paste0("geno",i),geno)
geno.all<- geno
seqpos<- round(as.numeric(colnames(geno.all))*NBASES)

#---------------------
# Calculate summary statistics
#---------------------
	 # JUST ONE SAMPLE
	# Combine to get genotypes
	gt <- Rcpp_make_single(geno.all)
	#NUMBER OF SUMSTATS
	#this is going tom be used to calculate r^2 dist see below
	testpos = seq(2.5,log10(NBASES),by=0.3) # to have into account the size of the fragment
	testpos = testpos[-length(testpos)] #remove last position
	countr2=length(testpos)

	nb_ss<-1+9+10+countr2+16

	sumstat <- rep(NA, nb_ss)
	
	#pi
	#sumstat[1] <- pairwise(geno.all) # pi NOT USED (real data is in 012 format)
	#count=2

	# Number of segregating sites (all pops
	sumstat[1] <- length(geno.all[1, ])
	count=2
	
	#pairwise euclidean distance for genotype matrix (0,1,2)
	#this is for: quantiles c(0.01,0.05,0.1,0.2,0.5,0.8,0.9,0.95,0.99)     #TODO # Summaries of pairwise distances for genotypes  7columns: mean, sd, quantiles(0.05,0.25,0.5,0.75,0.95)
	sumstat[count:(count+8)] <- pairwise.sing_compiled(gt)[-c(1,2)] 
	count=count+9
	
	#fspec calc (currently assuming ancestral info not used (i.e. wrapped) and more than 11 individuals in the sample - otherwise
	#the function will just return the wrapped fspec). TODO# Summarise site frequency spectrum: quantiles (0.01,0.05,0.25,0.5,0.75,0.95,0.99)
	
	sumstat[count:(count+9)] <- summarise_fspec_compiled(fspec_compiled(gt))[,2] # Summarise site frequency spectrum: 10
	count=count+10

	#r^2 calc
	#r2dist(NumericMatrix seqs, double varthresh,NumericVector pos,NumericVector cutpos, double cutthresh,double maxthresh){
	# seqs is the alignment
	# varthresh is the threshold variance. A cutoff of 0.1 corresponds to 0.05 < p < 0.95
	# pos contains the base positions of the snps
	# cutpos is the vector of positions at which to measure R^2 (in log10)
	# cutthresh is the tolerance as a proportion of cutpos within which R^2 is calculated - e.g. 0.1 means +/- 10% of target position
	#              this is measure on the natural (not log) scale
	# maxthresh is the maximum window size - because 10% on a scale of MB contains many more SNPs than on a scale of 100bp
	# E.g. use 1000

	#testpos = seq(2.5,6,by=0.5)
	
	sumstat[count:(count+countr2-1)] <- r2dist(gt,0.1,seqpos,testpos,0.1,1000)
	count=count+countr2

	#roh calc
	#0.1 in examp above is for the same varthresh in r2dist (doesn't need to be the same value; 0 probably OK; but might
	#	be worth having to protect against rare variants that are errors)
	#6 in examp above (log10 scale) should be set to be the comparable with the max of the roh for the target data
	#needs to be the same for all simulated data to compare like with like
	res3 = HH_fun2_compiled(gt,seqpos,0.1)
	#print(round(res3[,2],2))
	# res3[!is.finite(res3)] = -20.0 #I'm a bit worried whether this is OK #there may be no very long roh 
					#need to worry about what do with this. For target need to set max so that there are no infs
					#so then if sim data has infs/NA, we just remove as outside threshold.
	sumstat[count:(count+15)] <- res3[,2] # 16 quantilesc(0.005,0.01,0.02,0.05,0.1,0.2,0.5,0.75,0.9,0.975,0.99,0.9975,0.999,0.9995,0.99975,0.9999)


	


write.table(t(as.matrix(sumstat)), paste("summary_stats_",args[4],".txt",sep=""), sep="\t",col.names=F,row.names=F, append=T)


