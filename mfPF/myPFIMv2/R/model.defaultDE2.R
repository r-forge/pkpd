# Author: Martin Fink
# Date: March 2011
#
# model.defaultAE <- function()
# Get default model and options to run myPFIMv2
## Input
# none
## Return list
# PFIMmodel, Type, parModelName, parModel, parModelVar, parModelVarType, parObsName, parObsErr, parArms, time, mpOpt
##
model.defaultDE2 <- function() {

	default <- list()
	#formPK <- bolus_1cpt_Vk()[[1]]
	#formPD <- immed_lin_null(formPK)[[1]]
	
	default$PFIMmodel <- myDEfunc2
	default$Type <- 'DE'
	default$parArmsIC <- list(0,0)
	default$dosing <- list(cbind(c(0,24), c(30,30)), cbind(c(0,12,24,36), c(10,10,20,10)))
	
	default$parModelName <- c('V','k','Alin')
	default$parModel <- c(4.5,0.5,3)
	default$parModelVar <- c(0.2,0.4,0.15)
	default$parModelVarType <- 'exp'
	default$parObsName <- c('Conc','Effect')
	default$parObsErr <- list(c(0.3,0.5), c(0,0.3))
	default$parObsTimes <- list(list(c(0.5,1,2,4,8,11.9,23.9,35.9,47.9),c(0.5,1,2,8,11.9)),
							list(c(0.5,1,2,4,8,11.9,23.9,35.9,47.9),c(0.5,1,2,8,11.9,23.9,35.9,47.9)))
	default$parArmsName <- c('tst')
	default$parArms  <- list(c(5), c(2))
	default$ArmsName <- list('30 mg','10 mg')
	default$TimeName <- 'time (hr)'
	default$tRange <- c(0,48)
	default$mpOpt <- list()

	# Return value: example structure
	return(default)
}

myDEfunc2 <- function(t, y, parModel, parArms, mpOpt) {
#	form1<-paste("dose/V*(exp(-k*t))")
	V <- parModel[1]
	k <- parModel[2]
	Alin <- parModel[3]
	
	conc <- y[1]/V
	eff <- Alin*conc
	
	dy <- -k*y[1]

	return(list(c(dy),c(conc,eff)))
}
