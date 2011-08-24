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
model.defaultDEjac2 <- function() {

	default <- list()
	
	default$PFIMmodel <- myDEJfunc2
	default$Type <- 'DE'
	default$parArmsIC <- list(expression(c(10*70/0.15, RateT*V/CLT)))
	# Partial derivative of IC wrt parameters
	default$JacIC=list(expression(c(0,0, 0,RateT/CLT, 0,0, 0,-RateT*V/CLT^2, 0,0, 0,V/CLT)))
	# Partial derivative function for observations (if not equal to all states)
	default$JacObs=TRUE

	default$parModelName <- c("Kd","V","CLD","CLT","CLC","RateT")
	default$parModel <- c(0.2,6,0.2,1,0.4,40)
	default$parModelVar <- c(0.05,0.30,0.2,0.25,0.25,0.50)
	default$parModelVarType <- 'exp'
	default$parObsName <- c('TotDrug','TotTarget','Complex')
	default$parObsErr <- list(c(0,0.2), c(0,0.1), c(0.1,0.25))
	default$parObsTimes <- list(c(10/60/24,1/24,6/24,0.5, 2, 4, 7, 14, 21, 35, 49, 63, 77, 91), c(10/60/24, 7, 35, 77, 91), c(0.5, 2, 4, 7, 14, 21, 35, 49, 63, 77, 91))
	default$parArmsName <- c('dummy')
	default$parArms  <- list(c(0))
	default$ArmsName <- list('TMDD model')
	default$TimeName <- 'time (hr)'
	default$tRange <- c(0,91)
	default$mpOpt <- list()
	default$subjects <- c(100)

	# Return value: example structure
	return(default)
}

myDEJfunc2<-function(t, y, p, parArms, mpOpt){
	Kd  <- p[1]
	V   <- p[2]
	CLD <- p[3]
	CLT <- p[4]
	CLC <- p[5]
	RateT<-p[6]
	TD <- y[1]
	TT <- y[2]

	C =((Kd*V+TD+TT)-sqrt((Kd*V+TD+TT)^2 -4*TD*TT))/2
	D = TD-C        #free drug: total drug minus complex
	T = TT-C        #free target: total target minus complex

	dTD <-       - D*CLD/V - C*CLC/V	# Total Drug
	dTT <- RateT - T*CLT/V - C*CLC/V	# Total Target

	CFD = 0.15*D/V        #µg/mL total antibody is free plus complex
	CTT = 190 *TT/V       #ng/mL free plus complex for total target
	CFT = 190 *T/V        #ng/mL free target

	## Patrial derivatives w.r.t parameters
	# Partial derivatives of C wrt param
	GradC <- c((V - 0.5*(2*(V*(Kd*V + TD + TT))*((Kd*V + TD + TT)^2 - 4*TD*TT)^-0.5))/2,
			(Kd - 0.5*(2*(Kd*(Kd*V + TD + TT))*((Kd*V + TD + TT)^2 - 4*TD*TT)^-0.5))/2, 0,0,0,0)
	# Partial derivatives of TD and TT wrt parameters
	GradTD <- c(0, D*CLD/V^2 + C*CLC/V^2, -D/V, 0, -C/V, 0) + GradC*(CLD-CLC)/V
	GradTT <- c(0, T*CLT/V^2 + C*CLC/V^2, 0, -T/V, -C/V, 1) + GradC*(CLT-CLC)/V
	JacParam <- rbind(GradTD, GradTT)

	## Patrial derivatives w.r.t states
	# Partial derivatives of C wrt states
	GradCstate <- c((1 - 0.5 * ((2*(Kd*V + TD + TT) - 4*TT) * ((Kd* V + TD + TT)^2 - 4*TD*TT)^-0.5))/2,
			(1 - 0.5*((2*(Kd*V + TD + TT) - 4*TD)*((Kd*V + TD + TT)^2 - 4*TD*TT)^-0.5))/2)
	# Partial derivatives of TD and TT wrt states
	GradTDstate <- c(-CLD/V, 0) + GradCstate*(CLD-CLC)/V
	GradTTstate <- c(0, -CLT/V) + GradCstate*(CLT-CLC)/V
	JacState <- rbind(GradTDstate, GradTTstate)

	## Patrial derivatives of observations w.r.t. param & states
	# Partial derivatives wrt parameters
	JacObsParam <- rbind(0*GradC, 0*GradC, GradC)
	# Partial derivatives wrt states
	JacObsState <- rbind(c(1,0),c(0,1),GradCstate)
	
	return(list(c(dTD,dTT), c(TD,TT, C), JacState, JacParam, JacObsParam, JacObsState))
}
