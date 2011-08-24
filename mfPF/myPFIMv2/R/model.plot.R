# Author: Martin Fink
# Date: March 2011
#
# model.plot <- function(model, model.check=TRUE, plotSens=F)
# Solves a model (analytical or ODE) for given parameters and sampling time points
## Input
# model  .. model specification
# model.check .. run model.check as well
# plotSens .. plot not only solution, but also sensitivities
## Return list
# Solution matrix (as model.solve or model.sens)
##
model.plot <- function(model, model.check=TRUE, plotSens=FALSE, res=NULL, ...) {
	# Do model check if specified (= not opted out)
	if (model.check==TRUE) model.check(model)

	if (!is.null(res)) {
	# If solution matrix is given, then check some basics...
		dimRes <- dim(res)
		if (dimRes[1] != length(model$fullTime)) stop("Given results matrix does not conform with length of model$fullTime vector.")
		if (plotSens == TRUE & dimRes[3]!=model$nModpar+1) stop("Sensitivities have not been included in the results matrix (should be 3-dimensional).")
	} else {
	# Otherwise: Run either solve or sens depending on requested output
		if (plotSens==FALSE) {
			res <- model.solve(model, model$fullTime)
		} else {
			res <- model.sens(model, model$fullTime)
		}
	}

	## Output
	# Define some parameters for user-friendliness
	nModpar = model$nModpar
	nObs = model$nObs
	nArms = model$nArms

	# ToDo:
	# Safeguards for large nModpar & nArms
	# Include labels, etc.
	# Use .plot.res <- function below
	par(mfrow=c(nObs,nArms))
	
	# Sensitivities
	if (plotSens==TRUE) {
		# First dimension of res-array distinguishes solution and partial derivatives
		# Second dimention provides all obs first arm, then all obs second arm, ...
		# Third dimension is for the time
		for (iObs in seq(nObs)) {
			for (iArms in seq(nArms)) {
				.plot.res(1+seq(nModpar), iObs, iArms, res, model, ...)
			}
		}
		return(res)
	}
	
	# Solution
		# First dimension of res-array distinguishes solution and partial derivatives
		# Second dimention provides all obs first arm, then all obs second arm, ...
		# Third dimension is for the time
	for (iObs in seq(nObs)) {
		for (iArms in seq(nArms)) {
			.plot.res(1, iObs, iArms, res, model, ...)
		}
	}

	# Return value of model.solve or model.sens
	return(res)
}

.plot.res <- function(sModpar, iObs, iArms, res, model, col='black', col.line=col, col.symbol='red', pch=1, ...) {
	# Add group/arm title on top of page for each arm
	if (iObs == 1) title = model$ArmsName[iArms] else title = NULL
	# For plotting sensitivities, add that to label
	if (sModpar[1]>1) model$parObsName <- paste('Partial derivatives of', model$parObsName)

	# Select the data columns to plot
	colIdx <- iObs + (iArms-1)*model$nObs
	
	# Obtain limits for y-axis
	thisLog <- ifelse(sModpar[1]==1 & model$parObsLogY[iObs], "y", "")
	if (sModpar[1]==1 & model$parObsLogY[iObs]) {
	# If logarithmic, then min(yLim) >= parObsLogYMin (after removing zeroes)
		res[ res[, colIdx, sModpar]==0  , colIdx, sModpar] <- NA
		yLim = c(max(min(res[, colIdx, sModpar], na.rm=T),model$parObsLogYMin[iObs]), max(res[, colIdx, sModpar], na.rm=T))
	} else {
	# Default
		yLim = c(min(res[, colIdx, sModpar]), max(res[, colIdx, sModpar]))
	}
	
	# Plot sampling time points or several lines for sensitivities
	plot(model$fullTime, res[, colIdx, sModpar[1]], type='l', col=col.line, xlab=model$TimeName, ylab=model$parObsName[iObs], ylim=yLim, main=title, log=thisLog, ...)

	if (sModpar[1]==1) {
	# Add sampling time points
		tim <- model$timeMap$time[model$timeMap$iObs==iObs & model$timeMap[,2+iArms]]
		lines(tim, res[model$fullTime %in% tim, colIdx, sModpar[1]], type='p', col=col.symbol, pch=pch, ...)
	} else {
	# Add the remaining sensitivities
		for (iModpar in sModpar[2:length(sModpar)])
			lines(model$fullTime, res[, colIdx, iModpar], type='l', col=iModpar, ...)
	}
}
