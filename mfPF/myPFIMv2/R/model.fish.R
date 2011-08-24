# Author: Martin Fink
# Date: June 2011
#
# model.fish <- function(model, mysens=NULL)
# Derives partial derivatives of solution to the model parameters for a given model (analytical or ODE) at the specified time points
## Input
# model  .. model specification
# mysens .. sensitivities derived by model.sens (if not provided or inconsistent it is re-calculated)
## Return list
# 3-dim solution matrix (2D symmetric fisher information matrix)
# Third dimension: number of arms
# Size: nModpar x nModpar x Arms
##
model.fish <- function(model, mysens=NULL) {
	# If mysens not given - recalculate
	if (is.null(mysens)) mysens <- model.sens(model)

	# Define user-friendly parameters
	nModpar <- model$nModpar
	nObs <- model$nObs
	nArms <- ncol(model$parArms)

	res <- list()
	for (iArm in 1:nArms) {
		# Select observations times for this arm
		obsIdx <- model$timeMap$iObs

		# Select the sampling time points
		selectTime <- model$timeMap[, 2+iArm]
		
		# Sensitivity array dimensions
		#mysens <- array(dim=c(length(time), 1+nModpar, nObs*nArms))
		# Select all observations from one arm
		idx <- (1+(iArm-1)*nObs):(iArm*nObs)
		thisSens <- mysens[,idx,,drop=FALSE]

		# Reshape solution and sensitivities by aligning all observations and sampling times in one dimension
		mod <- as.vector(thisSens[,,1]) ## model solution
		sensi <- thisSens[,,(1:nModpar)+1] ## sensitivities nModpar x nObs x nTime
		dim(sensi) <- c(length(mod), nModpar)

		# Derive fisher information matrix for this arm
		temp <- model.fish.bloc(mod[selectTime], t(sensi[selectTime,]), obsIdx[selectTime], model)
		res[[iArm]] <- temp[[1]]
		debug <- temp
	}
	return(list(res,debug))
}
