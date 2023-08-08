#NEW sumstat file 12 Nov 2019

############################################
# nsam is number of sampled genomes
read.ms.output <- function (txt, nsam) 
{   
    h <- numeric()
    result <- list()
    gamlist <- list()
    positions <- list()
    times <- NaN
    probs <- NaN
    marker <- grep("segsites", txt)
    ndraws <- length(marker)
    stopifnot(length(marker) == ndraws)
    cat("ndraw =",ndraws," & nsam =",nsam,"\n")
	segsites <- sapply(strsplit(txt[marker], split = " "), function(vec) as.integer(vec[2]))
    for (draw in seq(along = marker)) {
        if (segsites[draw] > 0) {
            tpos <- strsplit(txt[marker[draw] + 1], split = " ")
            positions[[draw]] <- as.numeric(tpos[[1]][2:(segsites[draw] + 
                1)])
            haplotypes <- txt[(marker[draw] + 2):(marker[draw] + 
                2 + nsam[draw] - 1)]
            haplotypes <- strsplit(haplotypes, split = "")
            h <- sapply(haplotypes, function(el) c(as.integer(el)))
            if (segsites[draw] == 1) 
                h <- as.matrix(h)
            else h <- t(h)
        }
        else {
            h <- matrix(nrow = nsam[draw], ncol = 0)
            positions[[draw]] <- NA
        }
        gamlist[[draw]] <- h
        stopifnot(all(dim(h) == c(nsam[draw], segsites[draw])))
    }
    list(segsites = segsites, gametes = gamlist, probs = probs, 
        times = t(times), positions = positions, nsam = nsam, 
        nreps = ndraws)
}

#####################################################################

pairwise = function(m1)
#this is for 0,1 haplotype data (e.g. straight out of ms)
{
	d1 = dist(m1,method="man") #we use a manhattan distance because this corresponds 
								#the pairwise difference 
	mean(d1) #as long as we use a manhattan distance the mean is unaffected by phase 
	         #but not any higher moments, which are strongly affected.
}

make.single = function(a)
#this is to make 0,1,2 diploid data
{
	nhap = length(a[,1])
	nind = floor(nhap/2)
	if(nhap != 2*nind)stop("make.single: number of haplotypes is not even")
	nsnp = length(a[1,])
	h1 = a[seq(1,nhap,by=2),]
	h2 = a[seq(2,nhap,by=2),]
	s1 = h1 + h2
	s1
}

pairwise.sing = function(s,scale=T)
#this uses genotype (0,1,2) diploid data, as from make.single().
{
	if(scale){
		s = scale(s)
	}
	d1 = dist(s) #this is for the genotype (0,1,2) data. We are going to use this
	             #to get summary of distribution. We are happy with standard Euclidean
	             #distance for this
	a1 = mean(d1)
	a2 = sd(d1)
	a3 = quantile(d1,c(0.01,0.05,0.1,0.2,0.5,0.8,0.9,0.95,0.99),na.rm=T)
	as.numeric(c(a1,a2,a3))
	
}

fspec = function(s,wrap=T)
{
	nsnp = length(s[1,])
	nhap = length(s[,1])*2
	freqs = colSums(s)
	t1 = table(freqs)
	#return(t1)
	f1 = as.numeric(names(t1))
	c1 = as.numeric(t1)
	if(min(f1) == 0 || max(f1) == nhap){ #we need to remove the monomorphic ones
		if(min(f1) == 0){
			f1 = f1[-1]
			c1 = c1[-1]	
		}
		if(max(f1) == nhap){
			f1 = f1[-length(f1)]
			c1 = c1[-length(c1)]	
		}
		if(length(f1)!=length(c1))stop("fspec: f1,c1 have different lengths")
		if(length(f1) > nhap-1)stop("fspec: f1,c1 too long")
	}
	f2 = c(1:(nhap-1))
	c2 = rep(0,nhap-1)
	c2[f1] = c1
	if(wrap){
		c3 = c2 + rev(c2)
		c3[nhap/2] = c3[nhap/2]/2 #because this is the one it wraps around, so duplicated
		return(c3[1:(nhap/2)])
	}
	return(c2)
	
}

old_summarise_fspec = function(v)
#don't use
{
	ip = length(v)
	if(ip <= 10)stop("summarise_fspec: currently assumes spectrum has length greater than 10")
	q1 = ceiling(rep(ip,6)/c(ip,10,5,2.5,1.75,1.25))
	a1 = cumsum(v)/sum(v)
	a1[q1]
}

summarise_fspec = function(v)
{
	ip = length(v) #assuming it is wrapped, from diploids, it should be the same as the number of individuals
			# why? 0..2n, removing monomorphic gives 1..(2n - 1) equals 2n-1 - 1 + 1 = 2n-1. Always odd, 
			#so wrap around n
	cumdist = cumsum(v)
	cumdist = cumdist/cumdist[ip]
	if(ip <= 11) return(cbind(c(1:(ip-1)),cumdist[-ip]))
	seq1 = seq(0,ip,length = 12)
	seq1 = floor(seq1[c(2:11)])
	return(cbind(seq1,cumdist[seq1]))
}

getcovs = function(m1)
{
	if(dim(m1)[2] > 15000){
		print(paste("getcovs: size of cov matrix > 15000: ",dim(m1)[2]))
		return(c(NA,NA))
	}
	xx = cov(m1)
	utri = upper.tri(xx)
	c1 = abs(xx[utri])
	c(mean(c1),sd(c1))
}

testfun = function(ind1,inpos)
{
	isum = 0;
	for(j in 1:nrow(ind1)){
		wt = ind1[j,] == 1
		isum = isum + sum(wt)
		dvec = diff(inpos[wt])
		dvec = c(inpos[wt][1],dvec)
		if(j == 1)roh = dvec
		else roh = c(roh,dvec)
	}
	print(isum)
	roh
}

logit = function(x){log(x/(1-x))}

HH_fun = function(geno,pos,varthresh,ulim)
#ulim (log10 scale) should be set to be the max of the roh for the target data ??
#needs to be the same for all simulated data to compare like with like
{
	roh = rohcalc(geno,varthresh,pos)
	Fn = ecdf(roh)
	cutx = seq(0,ulim,length=16)
	#round this to integer numbers, because cdf jumps at integers (noticeable effect for small integers)
	cutx = log10(round(10^cutx))
	out1 = logit(1-Fn(10^cutx))
	return(cbind(cutx,out1))

}

HH_fun2 = function(geno,pos,varthresh)
#needs to be the same for all simulated data to compare like with like
{
	query = c(0.005,0.01,0.02,0.05,0.1,0.2,0.5,0.75,0.9,0.975,0.99,0.9975,0.999,0.9995,0.99975,0.9999)
	roh = rohcalc(geno,varthresh,pos) + 1 #do this to avoid -infs in the log()
	q1 = quantile(roh,query)
	return(cbind(logit(query),log(q1)))
}

make_gam_formula = function(pos,y.namelist,x.namelist)
#given variable position in y, and namelists (x,y), makes a formula for gam()
{
	extobj = paste(y.namelist[pos],"~ ",sep="")
	for(j in 1:length(x.namelist)){
	fobj = paste("s(",x.namelist[j],") + ",sep="") # you need the s() for GAM
	extobj = paste(extobj,fobj,sep="")
	}
	extobj = substr(extobj,1,nchar(extobj)-3) #removes the final " + "
	fmla = as.formula(extobj)
	return(fmla)
}


# Compiled versions
require(compiler)
enableJIT(3)

read.ms.output_compiled <- cmpfun(read.ms.output)
pairwise_compiled <- cmpfun(pairwise)
pairwise.sing_compiled <- cmpfun(pairwise.sing)
fspec_compiled <- cmpfun(fspec)
summarise_fspec_compiled <- cmpfun(summarise_fspec)
#getcovs_compiled <- cmpfun(getcovs)
HH_fun2_compiled<- cmpfun(HH_fun2)

#check if functions are compiled
is.compile <- function(func)
{
	# this function lets us know if a function has been byte-coded or not
	#If you have a better idea for how to do this - please let me know...
    if(class(func) != "function") stop("You need to enter a function")
    last_2_lines <- tail(capture.output(func),2)
    any(grepl("bytecode:", last_2_lines)) # returns TRUE if it finds the text "bytecode:" in any of the last two lines of the function"s print
}
	

