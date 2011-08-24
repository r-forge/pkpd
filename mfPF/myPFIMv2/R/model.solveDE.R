# Author: Martin Fink
# Date: March 2011
#
# model.solveDE <- function(parArms, time, mpOpt, model)
# Solves an ODE model for given parameters and time points
## Input
# model  .. model specification
# time      .. time points for output
## Return list
# Solution array (column vector for each observable)
# Size: length(time) x (nObs*nArms) x 1
##
model.solveDE <- function(model, time=model$time) {
	modelDE = model$PFIMmodel
	nModpar = model$nModpar
	nObs = model$nObs
	parArms = model$parArms
	nArms = model$nArms
	nArmspar = model$nArmspar
	parModel <- model$parModel

	rtol = .Machine$double.eps ^ 0.5 # default: 1e-6
	atol = .Machine$double.eps ^ 0.5 # default: 1e-6
	if (!is.null(model$Hmax)) hmax <- model$Hmax else hmax <- NULL
	
	# Set initial time to 0 if min(time) > 0
	setInitialTime = FALSE
	if (time[1] > 0) {
		time <- c(0, time)
		setInitialTime = TRUE
	}
	
	# Define general model parameters (for possible use in initial conditions)
	for (iModpar in seq(nModpar)) assign(model$parModelName[iModpar], parModel[iModpar])
	
	res <- array(dim=c(length(time), nObs*nArms, 1))
	for (iArms in seq(nArms)) {
		# Define arm/group parameters (for possible use in initial conditions)
		for (iArmpar in seq(nArmspar)) assign(model$parArmsName[iArmpar], parArms[iArmpar,iArms])
		IC <- eval(model$parArmsIC[[iArms]])
		obsIdx <- 1 + length(IC) + seq(nObs)
		# If dosing event entry exists, then add to lsoda
		if (!is.null(model$dosing)) {
			names(IC)[1] <- "v1"
			eventdat <- model$dosing[[iArms]]
		} else eventdat <- NULL
		res[, seq(nObs) + (iArms-1)*nObs, 1] <- lsoda(IC, time, modelDE, parModel, rtol=rtol, atol=atol, hmax=hmax, events = list(data = eventdat))[, obsIdx]
	}

	# If initial time was set to 0 (without having a measurement point there, we have to remove it again)
	if (TRUE == setInitialTime) {
		res <- res[-1, , , drop=F]
	}
	
	return(res)
}
