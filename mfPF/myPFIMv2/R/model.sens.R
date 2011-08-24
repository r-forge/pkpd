# Author: Martin Fink
# Date: March 2011
#
# model.sens <- function(model, parArms, time=NULL, useHess=FALSE)
# Derives partial derivatives of solution to the model parameters for a given model (analytical or ODE) at the specified time points
## Input
# model  .. model specification
# time      .. time points for output
# mpOpt  .. general option specification for myPFIMv2
## Return list
# 3-dim solution matrix (column vector for each observable)
# Third dimension: first layer for solution and then for each parModel the partial derivatives
# Size: (nModpar+1) x (nObs*nArms) x length(time)
##
model.sens <- function(model, time=model$time, useHess=FALSE) {

	# If algebraic equation branch out to analytical solution from AE solver
	if (model$Type == 'AE') {
		res = try(model.solveAE(model, time, calcSens=TRUE), silent=TRUE)
		# If analytical derivative was successful we can return - otherwise numerical approximation
		if (is.numeric(res)) return(res)
		warning('Could not derive analytical partial derivatives for the specified model. Now using numerical approximation.')
	}
	if (!is.null(model$JacIC)) {
		res = model.sensJac(model, time)
		return(res)
	}
	res = model.sensNum(model, time, useHess)
	return(res)
}
