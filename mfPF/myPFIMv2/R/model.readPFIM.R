# Author: Martin Fink
# Date: March 2011
#
# model.define <- function(PFIMmodel, Type, parModelName, parModel, parModelVar, parModelVarType, parObsName, parObsErr, parArmsName, TimeName, tInit, parArmsIC)
# Define a model (analytical or ODE) with inputs, outputs, ...
## Input
# stdin .. fileName for stdin PFIM file (including path)
## Return list
# model  .. model specification
##
model.readPFIM <- function(stdin, model_created=FALSE) {
	# Set PFIM program directory
	directory.program <- '../../tstPFIM3.2.1/Program/'
	examplePath <- dirname(stdin)
	assign('directory', examplePath, envir=.GlobalEnv)
	assign('dirsep', '/', envir=.GlobalEnv)
	
	# Read in default file
	source(stdin, local=T)
	source(paste(examplePath,'/',file.model,sep=''), local=T)
	mc_file_path <- paste(examplePath,'/model_created.r',sep='')
	if (model_created & file.exists(mc_file_path)) source(mc_file_path)

	## Current restrictions!
	if (covariate.model | covariate_occ.model | (n_occ>1))
		stop('Covariate model and IOV not implemented, yet.')
	if (run!='EVAL')
		stop('Only evaluation allowed at the moment.')
	if (option!=1)
		stop('Only block diagonal Fisher information matrix implemented.')

	## Collate given input parameters
	model = list()

	# AF or DE and respective equations
	model$Type <- modelform
	if (model$Type == 'AF') {
		model$Type <- 'AE'
		model$PFIMmodel <- form
	} else {
		model$PFIMmodel <- formED
	}
	
	# Model parameters
	model$parModelName <- parameters
	model$nModpar <- length(model$parModelName)
	model$parModel <- beta
	model$parModelVar <- diag(omega)
	model$parModelVarType <- Trand
	
	# Observations/measurements
	model$parObsName <- names.datay
	model$nObs <- nr
	if (model$nObs != length(model$parObsName))
		stop('names.datay should contain a label for all outputs/measurements/observations.')

	model$parObsErr <- list()
	model$parObsTimes <- list()
	for (idx in seq(model$nObs)) {
		model$parObsErr[[idx]] <- c(get(paste("sig.inter",LETTERS[idx],sep='')),get(paste("sig.slope",LETTERS[idx],sep='')))
		model$parObsTimes[[idx]] <- get(paste("prot",LETTERS[idx],sep=''))
	}
	model$parObsTimes$PFIM = TRUE

	if (exists('dose') & model$Type == 'AE') {
		model$parArms <- as.list(dose)
		model$parArmsName <- c('dose')
		model$ArmsName <- as.list(paste('Dose',dose))
	} else if (model$Type == 'DE') {
		nArms <- length(condinit)
		model$parArms <- rep(list(c(0)), nArms)
		model$parArmsName <- c('dummy')
		model$ArmsName <- as.list(paste('Output',seq(nArms)))
	} else {
		model$parArms <- list(c(0))
		model$parArmsName <- c('dummy')
		model$ArmsName <- list('Output')
	}
	model$nArms <- length(model$ArmsName)
	model$nArmspar <- length(model$parArmsName)

	model$TimeName <- names.datax[1]
	model$subjects <- subjects

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
			stop("Length of parameters in parArms list does not match number of arms/groups * number of parameters/arm.")
		model$parArms <- matrix(unlist(model$parArms), nrow=nParArms, ncol=nArms)
	}
	
	# For DE ensure that respective libraries are loaded
	if (model$Type == 'DE') {
		require(deSolve)
		if (nArms == 1)
			model$parArmsIC	<- list(condinit)
		else
			model$parArmsIC	<- as.list(condinit)
		model$tRange <- c(time.condinit, max(unlist(model$parObsTimes)))
		model$Hmax <- Hmax
	} else {
		model$tRange <- c(0, max(unlist(model$parObsTimes)))
	}
	
	if (model$parModelVarType %in% c(1,2)) {
		model$parModelVarType <- c('add','exp')[model$parModelVarType]
	}

	model <- model.define(model)
	model.check(model)
	
	# Return value: model structure
	return(model)	
}
