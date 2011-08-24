# Author: Martin Fink
# Date: March 2011
#
# model.sens <- function(model, parArms, time, mpOpt)
# Derives partial derivatives of solution to the model parameters for a given model (analytical or ODE) at the specified time points
## Input
# model  .. model specification
# time      .. time points for output
## Return list
# 3-dim solution matrix (column vector for each observable)
# Third dimension: first layer for solution and then for each parModel the partial derivatives
# Size: (nModpar+1) x (nObs*nArms) x length(time)
##
model.sensJac <- function(model, time=model$time) {

	# Define user-friendly parameters
	nModpar <- model$nModpar
	nObs <- model$nObs
	nArms <- model$nArms
	parModel <- model$parModel
	parArms <- model$parArms

	rtol = .Machine$double.eps ^ 0.5 # default: 1e-6
	atol = .Machine$double.eps ^ 0.5 # default: 1e-6
	if (!is.null(model$Hmax)) hmax <- model$Hmax else hmax <- NULL
	
	# Set initial time to 0 if min(time) > 0
	setInitialTime = FALSE
	if (time[1] > 0) {
		time <- c(0, time)
		setInitialTime = TRUE
	}
	
	# Define general model parameters (for use in initial conditions)
	for (iModpar in seq(nModpar)) assign(model$parModelName[iModpar], parModel[iModpar])
	
	res <- array(dim=c(length(time), nObs*nArms, 1+nModpar))
	for (iArms in seq(nArms)) {
		# Define arm/group parameters (for use in initial conditions)
		for (iArmpar in seq(model$nArmspar)) assign(model$parArmsName[iArmpar], parArms[iArmpar,iArms])

		# Derive initial conditions
		IC <- eval(model$parArmsIC[[iArms]])
		model$nIC <- length(IC)

		# Include ICs for partial derivatives
		if (length(JacIC<-eval(model$JacIC[[iArms]])) == 0)
			stop("Analytical sensitivities with ODES selected: Supposed to provide partial differentials, but missing model$JacIC, derivatives of initial conditions")
		IC <- c(IC, JacIC)

		# If dosing event entry exists, then add to lsoda
		if (!is.null(model$dosing)) {
			names(IC)[1] <- "v1"
			eventdat <- model$dosing[[iArms]]
		} else eventdat <- NULL

		# Solve ODE system (with analytical partial derivatives)
		model$parArmsNow <- model$parArms[[iArms]]
		tmp <- lsoda(IC, time, .model.solve.jacUser.interface, model, rtol=rtol, atol=atol, hmax=hmax, events = list(data = eventdat))[,-1]

		# If observations are not the state variables but derived variables
		if (dim(tmp)[2] != nObs*(1+nModpar)) {
			tmp <- tmp[, model$nIC*(1+nModpar)+seq(nObs*(1+nModpar))]
		}

		dim(tmp) <- c(length(time), nObs, 1+nModpar)
		res[, seq(nObs) + (iArms-1)*nObs,] <- tmp
	}
	# If initial time was set to 0 (without having a measurement point there, we have to remove it again)
	if (TRUE == setInitialTime) res <- res[-1,,,drop=F]

	return(res)
}

.model.solve.jacUser.interface <- function(tim, y, model) {
	
	tmp <- model$PFIMmodel(tim, y[seq(model$nIC)], model$parModel, model$parArmsNow)
	# Normal right hand side
	dy <- tmp[[1]]
	# Partial derivatives...
	# dxi/dpj -> d/dt(dxi/dpj) = d/dpj(dxi/dt) =
	#                d/dpj(fi) = Sum[@fi/@xk*dxk/dpj] + @fi/@pj
	JacState <- tmp[[3]]; JacParam <- tmp[[4]]
	pDer <- y[model$nIC+seq(model$nIC*model$nModpar)]; dim(pDer) <- c(model$nIC,model$nModpar)
	dpDer <- JacState %*% pDer + JacParam

	# If observations are not just the state vectors
	Obs <- NULL
	dObs <- NULL
	if (length(tmp)>4) {
		# Solution vector for observations
		Obs <- tmp[[2]]
		# Partial derivatives of observations
		dObs <- tmp[[6]] %*% pDer + tmp[[5]]
	}
	
	return(list(c(dy, dpDer), Obs, dObs))
}
