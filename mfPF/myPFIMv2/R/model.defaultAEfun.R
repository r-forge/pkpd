# Author: Martin Fink
# Date: March 2011
#
# model.defaultAE <- function()
# Get default model and options to run myPFIMv2
#
# Note the descrepancy between 't' in the model code and 'default$time' for denoting the time variable
## Input
# none
## Return list
# PFIMmodel, Type, parModelName, parModel, parModelVar, parModelVarType, parObsName, parObsErr, parArms, time, mpOpt
##
model.defaultAEfun <- function() {

	default <- list()
	default$PFIMmodel <- myAEfunction
	default$Type <- 'AE'

	default$parModelName <- c('V','k','Alin')
	default$parModel <- c(4.5,0.5,3)
	default$parModelVar <- c(0.2,0.4,0.15)
	default$parModelVarType <- 'exp'
	default$parObsName <- c('Conc','Effect')
	default$parObsErr <- list(c(0.3,0.5), c(0,0.3))
	default$parObsTimes <- list(list(c(0.5,1,2,4,8,12,24,48),c(0,4,8,12,24,48)),
								list(c(0.5,1,2,4,8,12,24),c(0,4,8,12,24)),
								list(c(0.5,1,4,12,23.9,47.9,71.9),c(0,4,12,23.9,47.9,71.9)))
	default$parArmsName <- c('dose','nDoses')
	default$parArms  <- list(c(100,30,5), c(1,1,3))
	default$ArmsName <- list('100 mg','30 mg','10 mg')
	default$TimeName <- 'Time (hr)'
	default$tRange <- c(0,72)
	default$mpOpt <- list()

	# Return value: example structure
	return(default)
}

myAEfunction <- function(tim, parModel, parArms, model) {
	V <- parModel[1]
	k <- parModel[2]
	Alin <- parModel[3]
	dose <- parArms[1]
	nDoses <- parArms[2]
	
	PK <- 0
	for (idx in 0:(nDoses-1)) {
		PK <- PK + (tim>=idx*24)*dose/V*(exp(-k*(tim-idx*24)))
	}
	PD <- Alin*PK

	return(cbind(PK,PD))
}
