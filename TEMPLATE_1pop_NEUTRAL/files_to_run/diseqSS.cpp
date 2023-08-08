#include <RcppArmadillo.h>
using namespace Rcpp;
//#include "Rcpp.h"
//[[Rcpp::depends(RcppArmadillo)]]


//[[Rcpp::export]]
NumericVector
r2dist(NumericMatrix seqs, double varthresh,NumericVector pos,NumericVector cutpos, double cutthresh,double maxthresh){
// seqs is the alignment
// varthresh is the threshold variance. A cutoff of 0.1 corresponds to 0.05 < p < 0.95
// pos contains the base positions of the snps
// cutpos is the vector of positions at which to measure R^2 (in log10)
// cutthresh is the tolerance as a proportion of cutpos within which R^2 is calculated - e.g. 0.1 means +/- 10% of target position
//              this is measure on the natural (not log) scale
// maxthresh is the maximum window size - because 10% on a scale of MB contains many more SNPs than on a scale of 100bp
// E.g. use 1000

	int nind = seqs.nrow(); //number of individuals
	int slen = seqs.ncol(); //length of sequence
	//Rcout << nind << "\n"; //use this to print out data during run
	//Rcout << slen << "\n"; //use this to print out data during run
	
	if(slen != pos.size())stop("length of position vector not same as num cols in sequence alignment");
	int ncut = cutpos.size();
	//Rcout << ncut << "\n"; //use this to print out data during run
	NumericVector mincut(ncut),maxcut(ncut);
	for(int i=0;i<ncut;++i){
		mincut[i] = pow(10,cutpos[i])*(1-cutthresh);
		maxcut[i] = pow(10,cutpos[i])*(1+cutthresh);
		if(maxcut[i] - mincut[i] > maxthresh){
			mincut[i] = pow(10,cutpos[i]) - maxthresh/2;
			maxcut[i] = pow(10,cutpos[i]) + maxthresh/2;
		}
	}

	IntegerVector icount(ncut);
	NumericVector r2sum(ncut);
	NumericVector varvec(slen);
	for(int i = 0;i < slen; i++)varvec[i] = var(seqs(_, i));
	for(int k = 0;k < ncut;++k){icount[k] = 0;r2sum[k] = 0.0;}
	//slen = 2000;
	for(int i = 0; i < slen;i++){
		if(varvec[i] < varthresh)continue;
		for(int j = i+1;j < slen-1;j++){
			if(varvec[j] < varthresh)continue;
			double dist = pos[j] - pos[i];
			int k;
			for(k = 0;k < ncut; ++k){
				if(dist > mincut[k] && dist < maxcut[k])break;
			}
			if(k == ncut){//not in useful distance, so discard
				continue;
			}
			NumericVector zzcol1 = seqs(_,i);
			NumericVector zzcol2 = seqs(_,j);
			arma::vec adash = as<arma::vec>(zzcol1);
			arma::vec bdash = as<arma::vec>(zzcol2);
			double temp1 =  arma::as_scalar(arma::cor(adash,bdash));
			//double temp1 = mycor(as<NumericVector>(zzcol1),as<NumericVector>(zzcol2),varthresh);
			if(NumericVector::is_na(temp1))continue;
			icount[k] = icount[k] + 1;
			r2sum[k] = r2sum[k] + pow(temp1,2.0);
			//Rcout << k << " " << dist << " " << icount[k] << " " << r2sum[k] << "\n";
			// to make it more efficient add a filter in this func to weed out cols < varthresh
		}
	}
	for(int k = 0;k < ncut;++k){
		if(icount[k] == 0){r2sum[k] = NA_REAL;}
		else r2sum[k] = r2sum[k]/icount[k];
		//Rcout << icount[k] << " ";	
	}
//	Rcout << "\n";
//	Rcout << sum(icount) << "\n";	
	r2sum = log(r2sum/(1-r2sum));
	return r2sum;
			
}


//[[Rcpp::export]]
NumericVector
rohcalc(NumericMatrix geno, double varthresh,NumericVector pos){

	int nind = geno.nrow(); //number of individuals
	int nsnp = geno.ncol(); //length of sequence

	NumericVector varvec(nsnp);
	for(int i = 0;i < nsnp; i++)varvec[i] = var(geno(_, i));

	LogicalMatrix wtmat(nind,nsnp);
	for(int j=0;j<nind;++j){
		for(int i=0;i<nsnp;++i){
			wtmat(j,i) = geno(j,i) == 1 && varvec[i] >= varthresh;
		}
	}
	int nitem = sum(wtmat);
	//Rcout << nitem << "\n";
// so now we need to think about how many segments this equals
// does it start at 0; I think it does, and R code assumes this at the moment.
// I think the vector needs to be nitem long;
	NumericVector outvec(nitem);
	int icount = 0;
	for(int j=0;j<nind;++j){
		double lastpos = 0;
		for(int i=0;i<nsnp;++i){
			if(wtmat(j,i)){
				outvec[icount++] = pos[i] - lastpos;
				lastpos = pos[i];
			}
		}
	}
	//Rcout << icount << "\n";
	return outvec;
}


//[[Rcpp::export]]
NumericMatrix Rcpp_make_single(NumericMatrix seqs){

	int nhap = seqs.nrow(); //number of individuals
	int nsnp = seqs.ncol(); //length of sequence
	int nind = floor(nhap/2);
	if(nhap != 2*nind)stop("make.single: number of haplotypes is not even");

	NumericMatrix indiv(nind,nsnp);
	for(int j=0;j<nind;++j){
		int j1 = 2*j;
		int j2 = j1+1;
		indiv(j,_) = seqs(j1,_) + seqs(j2,_);
	}
	return indiv;
}

