# Author: Martin Fink
# Date: March 2011
#
# model.solveAE <- function(parArms, time, model, mpOpt, calcSens=F)
# Solves an analytical model for given parameters and time points
# and/or derives the partial derivatives w.r.t. the parameters
## Input
# model  .. model specification
# time      .. time points for output
# calcSens  .. calculate patrial derivatives (sensitivities)
## Return list
# Solution array (column vector for each observable)
# Size: 1 x (nObs*nArms) x length(time) or (1+nModpar) x (nObs*nArms) x length(time)
##
model.solveAE <- function(model, tim=model$time, calcSens=F) {
	nModpar = model$nModpar
	parModel = model$parModel
	nObs = model$nObs
	parArms = model$parArms
	nArms = model$nArms
	nArmspar = model$nArmspar

	if (calcSens==F) {
		res <- array(dim=c(length(tim), nObs*nArms, 1))
	} else {
		res <- array(dim=c(length(tim), nObs*nArms, 1+nModpar))
	}

	# Assign model parameter values
	for (iModpar in seq(nModpar)) assign(model$parModelName[iModpar], parModel[iModpar])

	# For PFIM-like expression as model
	if (is.expression(model$PFIMmodel)) {
		t <- tim
		modelEQ <- model$PFIMmodel
		for (iArms in seq(nArms)) {
			for (iArmpar in seq(nArmspar)) assign(model$parArmsName[iArmpar], parArms[iArmpar,iArms])
			for (iObs in seq(nObs)) {
				res[, iObs + (iArms-1)*nObs, 1] <- eval(modelEQ[iObs])
			}
		}

		if (calcSens==TRUE) {
			for (iArms in seq(nArms)) {
				for (iArmpar in seq(nArmspar)) assign(model$parArmsName[iArmpar], parArms[iArmpar,iArms])
				for (iModpar in seq(nModpar)) {
					for (iObs in seq(nObs)) {
						modelEQ <- D(model$PFIMmodel[iObs], model$parModelName[iModpar])
						res[, iObs + (iArms-1)*nObs, 1+iModpar] <- eval(modelEQ)
					}
				}
			}
		}
		return(res)
	}
	
	# For model <- function(tim, parModel, parArms, opt)
	if (is.function(model$PFIMmodel)) {
		# If sensitivities should be calculated then do this by numerical derivatives
		if (calcSens==TRUE) return(model.sensNum(model, tim))
		# This provides the solution of the system by function call
		for (iArms in seq(nArms)) {
			tmp = try(withCallingHandlers(model$PFIMmodel(tim, parModel, parArms[,iArms], model), warning = stop()), silent=TRUE)
			if (is.numeric(tmp) & length(tmp)==length(tim)*nObs) {
				# If the function is 'vectorized' it immediately prints out all observations at all times
				res[,(iArms-1)*nObs + seq(nObs),1] <- tmp
			} else {
				# If the function is not as flexible we have to cycle through all time points
				res[,(iArms-1)*nObs + seq(nObs),1] <- t(sapply(1:length(tim), function(tidx, parModel, parA, model) model$PFIMmodel(tim[tidx], parModel, parA, model), parModel=parModel, parA=parArms[,iArms], model=model))
			}
		}
		return(res)
	}
	
	warning("Could not determine model type (expression, function, ODE,...) properly")
	return(NULL)
}
