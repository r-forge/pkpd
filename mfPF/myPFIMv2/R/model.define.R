# Author: Martin Fink
# Date: March 2011
#
# model.define <- function(PFIMmodel, Type, parModelName, parModel, parModelVar, parModelVarType, parObsName, parObsErr, parArmsName, TimeName, tInit, parArmsIC)
# Define a model (analytical or ODE) with inputs, outputs, ...
## Input
# PFIMmodel .. Model definition derived from PFIM models
#	- Algebraic models given as 'expression' using the parameter names as specified in parModelName
#	- ODE models given as function(parModel, parArms, time, mpOpt)
# Type   .. AE or DE
# parModelName .. model parameter names
# parModel .. model parameter values
# parModelVar .. model parameter variances
# parModelVarType .. exponential or additive variance model
# parObsName .. observation/response names
# parObsErr .. error amount and structure related to observations
# parObsTimes .. possible sampling times per observation
# parArms .. arm/group-specific parameters
# parArmsName .. names of arm/group-specific parameters
# TimeName .. label for time axis
# parArmsIC .. initial conditions (for ODEs only)
## Return list
# model  .. model specification
##
model.define <- function(PFIMmodel, Type, parModelName, parModel, parModelVar, parModelVarType, parObsName, parObsErr, parObsTimes, parArms, parArmsName, ArmsName, TimeName, tRange, parArmsIC) {

	if (is.list(PFIMmodel)) {
		model <- PFIMmodel
	} else {
		# Collate given input parameters
		model = list()
		model$PFIMmodel <- PFIMmodel
		model$Type <- Type
		model$parModelName <- parModelName
		model$parModel <- parModel
		model$parModelVar <- parModelVar
		model$parModelVarType <- parModelVarType
		model$parObsName <- parObsName
		model$parObsErr <- parObsErr
		model$parObsTimes <- parObsTimes
		model$parArms <- parArms
		model$parArmsName <- parArmsName
		model$ArmsName <- ArmsName
		model$TimeName <- TimeName
		model$tRange <- tRange
		if (model$Type == 'DE')	model$parArmsIC	<- parArmsIC
	}

	# Number of parameters and observations
	model$nModpar <- length(model$parModelName)
	model$nObs <- length(model$parObsName)
	model$nArmspar <- length(model$parArmsName)
	model$nArms <- length(model$ArmsName)

	# Add info for graphical logarithmic output
	if (is.null(model$parObsLogY)) model$parObsLogY <- c(TRUE, rep(FALSE, model$nObs-1)) # Assume that first obs most often PK
	if (is.null(model$parObsLogYMin)) model$parObsLogYMin <- rep(1e-4, model$nObs)
	
	# Time components (time, observations, arms)
	if (is.null(model$dosing)) {
		model$time <- sort(unique(unlist(model$parObsTimes)))
	} else {
		tmp <- c() # add dosing time points to time vector for solution (but not observations)
		for (iArms in seq(model$nArms)) tmp <- c(tmp, model$dosing[[iArms]][,1])
		model$time <- sort(unique(c(unlist(model$parObsTimes), tmp)))
	}
	model$nTime <- length(model$time)
	model$fullTime <- sort(unique(c(seq(min(c(model$tRange[1], model$time[1])), max(model$tRange[2], model$time[2]), length.out=100), model$time)))
	tmp <- data.frame(time=rep(model$time, times=model$nObs), iObs=rep(1:model$nObs, each=model$nTime))
	if (!is.null(model$parObsTimes$PFIM)) { # PFIM - sorted by observations
		for (iObs in 1:(length(model$parObsTimes)-1)) {
			for (iArm in 1:length(model$parObsTimes[[iObs]])) {
				tmp[tmp$iObs==iObs, 2+iArm] <- (model$time %in% model$parObsTimes[[iObs]][[iArm]])
			}
		}
		if (model$nArms>1 & ncol(tmp)==3) {
		# If for all observations only time points for first arm are given - extend/copy to all arms
			tmp[,2+2:model$nArms] <- tmp[,3]
		}
	} else { # my structure - sorted by arm
		if (length(model$parObsTimes) < model$nArms)
			stop("Please provide observation times for each arm in model$parObsTimes.")
		for (iArm in 1:length(model$parObsTimes)) {
			if (is.list(model$parObsTimes[[iArm]])) {
				for (iObs in 1:length(model$parObsTimes[[iArm]])) {
					tmp[tmp$iObs==iObs, 2+iArm] <- (model$time %in% model$parObsTimes[[iArm]][[iObs]])
				}
			} else {
				# If for an arms only one set of time points are given - use/copy for all observations
				tmp[, 2+iArm] <- (model$time %in% model$parObsTimes[[iArm]])
			}
		}
	}
	model$timeMap <- tmp
	
	# Create matrices out of list objects (for being user friendly and avoid confusion!)
	if (is.list(model$parObsErr)) {
		if ( (nObsErr<-length(unlist(model$parObsErr))) %% 2 == 1)
			stop("Length of observation random error parameters must have 2 parameters per observation as list(c(add,prop),c(add,prop)).")
		model$parObsErr <- matrix(unlist(model$parObsErr), nrow=2, ncol=nObsErr/2)
	}
	if (is.list(model$parArms)) {
		nParArms <- length(model$parArmsName)
		nArms <- length(model$ArmsName)
		if (length(unlist(model$parArms)) != nParArms*nArms)
			stop("Length of parameters in parArms list does not match number of arms (groups) * number of parameters/arm.")
		model$parArms <- matrix(unlist(model$parArms), nrow=nParArms, ncol=nArms, byrow=TRUE)
	}
	
	# For DE ensure that respective libraries are loaded
	if (model$Type == 'DE')	{
		require(deSolve)
		if (!is.null(model$dosing))
			if (!is.data.frame(model$dosing)) {
			# If dosing is data.frame we assume that it is done according to lsoda needs...
			# Otherwise we implement all doses into the first compartment
				eventdat <- list()
				for (iArms in seq(nArms)) {
					eventdat[[iArms]] <- data.frame(var = "v1",
							time = model$dosing[[iArms]][,1],
							value = model$dosing[[iArms]][,2],
							method = "add")
				}
				model$dosing <- eventdat
			}
	}
	
	# PFIM backward compatibility
	if (model$parModelVarType %in% c(1,2)) {
		model$parModelVarType <- c('add','exp')[model$parModelVarType]
		warning("Random effect model should be either 'add' or 'exp'. ('1' and '2' currently allowed for backward compatibility)")
	}

	model.check(model)
	
	# Return value: model structure
	return(model)	
}
