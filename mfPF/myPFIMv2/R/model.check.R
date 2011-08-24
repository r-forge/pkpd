# Author: Martin Fink
# Date: March 2011
#
# model.check <- function(model, parArms, mpOpt)
# Checks model input and output consistency
# Runs the model file once for getting number of derivatives, etc.
## Input
# parArms  .. parameters for the different arms (incl. dosing, initial conditions,...)
# model  .. model specification
# mpOpt  .. general option specification for myPFIMv2
## Checks
# Run model once
# Check input consistency
# Check output consistency with specifications
## Output
# Produces error and warnings while running - no return value
##
model.check <- function(model) {

	# Model parameter specifications
	if (length(model$parModelName)!=length(model$parModel))
		stop("Length of parameter names of model and length of parameter values must coincide.")

	if (length(model$parModelName)!=length(model$parModelVar))
		stop("Length of parameter names of model and length of parameter variances must coincide.")

	if (any(model$parModel==0))
		stop("Model parameters can not be equal to zero. Please, reqrite your model accordingly.")

	if (!model$parModelVarType %in% c('exp','add'))
		stop("Random effect model must be either 'exp' or 'add'.")
		
	if (any(model$parModelVar<0))
		stop("Variances of model parameters must be positive.")

	# Observation specifications
	if (length(model$parObsName)!=ncol(model$parObsErr))
		stop("Length of observation names of model and number of columns of observation error specifications must coincide.")
	if (2!=nrow(model$parObsErr))
		stop("Need to specify additive and proportional error per observation in parObsErr.")
	if (length(model$parObsLogY)!=model$nObs) 
		stop("Length of choice of logarithmic output (parObsLogY) must be equal to number of observation types/names.")
	if (length(model$parObsLogYMin)!=model$nObs)
		stop("Length of minimum value for logarithmic output (parObsLogYMin) must be equal to number of observation types/names.")

	# Arms/groups specifications
	if (length(model$parArmsName)!=nrow(model$parArms))
		stop("Length of (arms/groups) parameter names and number of rows of their parameter matrix must coincide.")
	if (length(model$ArmsName)!=ncol(model$parArms))
		stop("Length of arms/groups names and number of columns of their parameter matrix must coincide.")
	
	# General model type, structure, etc.
	if (!model$Type %in% c('DE','AE'))
		stop("Model type must be either 'DE' or 'AE'.")

	if (model$Type == 'AE') {
		# Check length of observations with number of equations
		if (is.expression(model$PFIMmodel)) {
			if (model$nObs != length(model$PFIMmodel))
				stop("Number of observation names does not coincide with number of outputs in the algebraic model.")
		} else {
			if (model$nObs != length(model$PFIMmodel(0,model$parModel,model$parArms[,1],model)))
				stop("Number of observation names does not coincide with number of outputs in the algebraic model.")
		}
	} else { # 'DE'
		# Do something
		if (!is.null(model$dosing)) {
			if (model$nArms != length(model$dosing))
				stop("Number of entries in model$dosing must coincide with number of arms specified.")
			if (!is.data.frame(model$dosing[[1]]))
				for (iArms in seq(model$nArms)) {
					if (ncol(model$dosing[[iArms]]) != 2)
						stop("Columns of each entry in model$dosing must be for time and respective dose, i.e., 2 columns.")
				}
		}
		# Check validity of ICs
		#for (i in 1:p) {assign(parameters[i],theta[i]) }
		#cond<-eval(condinit[kit])	
	}
	
	# ToDo:
	# Check if all parameters exist in model
	# Run model once to check output structure
	return(TRUE)
}
