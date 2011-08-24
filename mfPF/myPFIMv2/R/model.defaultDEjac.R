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
model.defaultDEjac <- function() {

	default <- list()
	
	default$PFIMmodel <- myDEJfunc
	default$Type <- 'DE'
	default$parArmsIC <- list(expression(c(0,Rin/kout)))
	#default$parArmsIC <- list(expression(c(100/V,Rin/kout)))
	# Partial derivative of IC wrt parameters - indicates that other partial derivatives are returned by righthand side function as well (3rd & 4th return arguments)
	default$JacIC=list(expression(c(0,0, 0,0, 0,0, 0,1/kout, 0,-Rin/kout^2, 0,0)))

	default$parModelName <- c("V","Vm","km","Rin","kout","C50")
	default$parModel <- c(12,0.1,0.5,6.4,1.2,1)
	default$parModelVar <- c(0.25,0.25,0,0.3,0.25,0)
	default$parModelVarType <- 'exp'
	default$parObsName <- c('Conc','Effect')
	default$parObsErr <- list(c(0,0.2), c(3.8,0))
	default$parObsTimes <- list(c(0.5, 1, 2, 19, 38, 61, 160),c(0, 0.7, 1.5, 23, 12, 44, 144))
	default$parArmsName <- c('dummy')
	default$parArms  <- list(c(0))
	default$ArmsName <- list('PKPD model')
	default$TimeName <- 'time (hr)'
	default$tRange <- c(0,160)
	default$mpOpt <- list()
	default$subjects <- c(100)

	# Return value: example structure
	return(default)
}

myDEJfunc <- function(tim, y, parModel, parArms, mpOpt) {
	V<-parModel[1]
	Vm<-parModel[2]
	km<-parModel[3]
	Rin<-parModel[4]
	kout<-parModel[5]
	C50<-parModel[6]
	pk<-y[1]
	pd<-y[2]
	
	## PK
	dpk1<- - Vm*pk/(km*V + pk)
	# Partial derivative of PK wrt parameters
	GradPK <- c(Vm*pk/(km*V + pk)^2*km,
			- pk/(km*V + pk),
			+ Vm*pk/(km*V + pk)^2*V,
			0,0,0)
	# Partial derivative of PK wrt states
	GradPKstate <- c(-Vm/(km*V + pk) + Vm*pk/(km*V + pk)^2, 0)

	if(tim<=1){
		dpk1 <- dpk1 + 100/V
		GradPK[1] <- GradPK[1] -100/V^2
	}

	## PD
	dpd1 <- Rin*(1 - pk/(pk + C50)) - kout*pd
	# Partial derivative of PD wrt parameters
	GradPD <- c(0,0,0,
			(1 - pk/(pk + C50)),
			-pd,
			Rin*pk/(pk + C50)^2)
	# Partial derivative of PD wrt states
	GradPDstate <- c(-Rin/(pk + C50) + Rin*pk/(pk + C50)^2, -kout)

	## Partial derivatives of right hand side to states and parameters...
	Jac <- rbind(GradPK, GradPD)
	JacState <- rbind(GradPKstate, GradPDstate)				
	
	return(list(c(dpk1,dpd1),c(pk,pd),JacState,Jac))
}
