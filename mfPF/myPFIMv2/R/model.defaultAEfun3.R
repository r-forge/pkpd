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
model.defaultAEfun3 <- function() {

	default <- list()
	default$PFIMmodel <- myAEfunction3
	default$Type <- 'AE'

	default$parModelName <- c('V','k','keff','Emax','EC50')
	default$parModel <- c(4.5,0.15,3,2.3,3)
	default$parModelVar <- c(0.2,0.4,0.15,0.4,0.05)
	default$parModelVarType <- 'exp'
	default$parObsName <- c('Conc','Effect')
	default$parObsErr <- list(c(0.3,0.5), c(0.2,0.7))
	# Sampling times per arm - and within arm per observation/measurement
	default$parObsTimes <- list(c(0.5,1,2,4,8,12,24,48))
	default$parArmsName <- c('dose')
	default$parArms  <- list(c(200))
	default$ArmsName <- list('PKPD model')
	default$TimeName <- 'Time (hr)'
	default$tRange <- c(0,48)
	default$mpOpt <- list()

	# Return value: example structure
	return(default)
}

myAEfunction3 <- function(tim, parModel, parArms, model) {
	V <- parModel[1]
	k <- parModel[2]
	keff <- parModel[3]
	Emax <- parModel[4]
	EC50 <- parModel[5]
	dose <- parArms[1]
	
	if (tim <= 7) { # Appearance of ADAs at t=7
		PK = dose/V*(exp(-k*tim))
		temp = dose/V*k/(k-keff)*(exp(-keff*tim)-exp(-k*tim))
	} else {
		PK = 0
		temp = 0
	}
	Eff <- Emax*temp/(EC50 + temp)
	PD <- 1+Eff

	return(cbind(PK,PD))
}
