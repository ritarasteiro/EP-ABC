phicalc = function(r,Q)
#this is the log partition function. NB there are a number of wrong representations of this in the web/literature, but I think this is correct
#(see explanations in Evernote)
{
	muv = as.numeric(chol2inv(Q)%*%r) #possibly a bit inefficient - could be given as argument
	val = 0.5*(-as.numeric(determinant(Q/(2*pi))$mod) + muv %*% Q %*% muv)
	val

}

epanvol = function(dst,bandwidth,dimens) # This gives normalised k.d. weight for each point (remember to put in bandwidth^dimens term)
{
	x = dst^2/bandwidth^2
	svol = volsphere(dimens)
	val = 0.5/svol*(dimens+2)*(1-x)
	val = val/bandwidth^dimens
	val
}

volsphere = function(dimens)
#volume of sphere with radius 1 (need to multiply by radius^dimens for other radius values)
{
	val = pi^(dimens/2)/gamma(dimens/2 + 1)
	val
}
