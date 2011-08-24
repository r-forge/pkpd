# Author: Martin Fink
# Date: June 2011
#
# model.out.det_rse <- function(fish, model, output=FALSE)
# Derives determinant, criterion, se and r.s.e from given fisher information matrix
## Input
# model  .. model specification
# fish .. fisher information matrix
# output .. write out results TRUE/FALSE
## Return list
# Determinant
# Criterion
# Output data.frame with StdErr and RSE
##
model.out.det_rse <- function(model, fish, output=FALSE) {
	# Prepare names, etc. for output - remove variances=0 and resid err=0
	out <- data.frame(Name=c(model$parModelName, model$parModelName, paste(rep(model$parObsName,each=2),c('add','prop'))),
		Value=c(model$parModel, model$parModelVar, as.vector(model$parObsErr)), StdErr=0, RSE=0,
		parType=c(rep(c(1,2),each=model$nModpar), rep(3,length(model$parObsErr))) )
	names(out)[c(1,4)] <- c(' ','RSE(%)')
	 
	# Remove rows with variance of error == 0
	out <- out[out$Value !=0,]

	# Fisher information matrix: determinant and criterion
    detFish <- det(fish)
	crit <- (detFish)^(1/dim(fish)[1]) 

	# Derive standard errors and relative standard errors
	#invFish <- try(solve(as.matrix(fish)))
	#if (is.character(invFish)) {
	#	invFish <- myPseudoinverse(as.matrix(fish))
	#	warning("Fisher information matrix singular. Pseudoinverse used.")
	#}
	#out[,3] <- sqrt(diag(invFish))

	# Derive diagonal elements of inverse of symmetric matrix
	#	using eigendecomposition (more stable than solve)
	mEigen <- eigen(fish, symm=T);
	out[,3] <- sqrt((mEigen$vectors*mEigen$vectors) %*% (1/mEigen$values))

	out[,4] <- out[,3]/out[,2]*100
	
	res <- list(detFish, crit, out)
	if (output==FALSE) return(res)

	ncat('Population fisher information matrix')
	print(fish)
	ncat('Determinant and Criterion')
	print(detFish)
	print(crit)

	ncat('EXPECTED STANDARD ERRORS')
	ncat('Fixed effect parameters')
	print(out[out$parType==1,1:4], digits=3, row.names=F)

	ncat('Random effect parameters (IIV/BSV)')
	print(out[out$parType==2,1:4], digits=3, row.names=F)

	ncat('Residual error')
	print(out[out$parType==3,1:4], digits=3, row.names=F)	

	return(res)
}

ncat <- function(string) {
	cat('\n')
	cat(string)
	cat('\n')
	cat(rep('=',nchar(string)/2+1))
	cat('\n')
}
