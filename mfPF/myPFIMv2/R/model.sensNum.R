# Author: Martin Fink
# Date: March 2011
#
# model.sensNum <- function(model, parArms, time=NULL, useHess=FALSE)
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
model.sensNum <- function(model, time=model$time, useHess=FALSE) {

	# Define user-friendly parameters
	nModpar <- model$nModpar
	nObs <- model$nObs
	nArms <- model$nArms
	parModel <- model$parModel
	parArms <- model$parArms
	
	# Options for jacobian
	methods.args=list(r=4,eps=1e-4,d=1e-2,zero.tol=sqrt(.Machine$double.eps))
	# If higher accuracy of Jacobian is wanted, increase from r=4 to r=6, eps=1e-4

	res <- array(dim=c(length(time), nObs*nArms, 1+nModpar))
	res[,,1] <- model.solve(model, time)
	
	# Two different ways to derive the partial derivatives:
	# jacobian - library(numDeriv)
	# fdHess - library(nlme) - used by PFIM
	if (!useHess) { # jacobian
		require(numDeriv)
		df1<-jacobian(.model.solve.jac.interface, model$parModel, model=model, time=time, method.args=methods.args)
		for (iModpar in seq(nModpar)) {
			tmp <- df1[,iModpar]
			dim(tmp) <- c(length(time),model$nObs*nArms)
			res[,,1+iModpar] <- tmp
		}
	} else { # fdHess
		require(nlme)
		for (it in seq(time)) {
			for (iObs in seq(nObs)) {
				for (iArms in seq(nArms)) {
					tmp<-fdHess(parModel,.model.solve.fdHess.interface,model=model,t=time[it],iObs=iObs,iArms=iArms)$gradient
					res[it,iObs+(iArms-1)*nObs,1+seq(nModpar)] <- tmp
				}
			}
		}
	}
	return(res)
}

.model.solve.jac.interface <- function(parModel, model, time) {
	model$parModel <- parModel
	res <- model.solve(model, time)
	# Create one-dimensional vector for calculating Jacobian
	dim(res) <- length(time)*model$nObs*model$nArms
	return(as.vector(res))
}

.model.solve.fdHess.interface <- function(parModel, model, t, iObs, iArms) {
	model$parModel <- parModel
	res <- model.solve(model, c(model$tRange[1], t))
	return(res[2,iObs+(iArms-1)*model$nObs,1])
}
