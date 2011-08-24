# Author: Martin Fink
# Date: March 2011
#
# model.solveAE <- function(model, time=NULL)
# Solves a model (analytical or ODE) for given parameters and time points
## Input
# parModel .. basic model parameters
# model  .. model specification
# parArms  .. list of parameters for the different arms (incl. dosing, initial conditions,...)
# time      .. time points for output
# mpOpt  .. general option specification for myPFIMv2
## Return list
# Solution array (column vector for each observable)
# Size: 1 x (nObs*nArms) x length(time)
##
model.solve <- function(model, time=NULL) {
	
	if (model$Type == 'AE') {
		return(model.solveAE(model, time))
	} else {
		return(model.solveDE(model, time))
	}
}
