# Author: Martin Fink
# Date: March 2011
#
# model.defaultAE <- function()
# Get default model and options to run myPFIMv2
#
# Note the descrepancy between 't' in the model code and 'default$time' for denoting the time variable
## Input
# none
## Return list
# PFIMmodel, Type, parModelName, parModel, parModelVar, parModelVarType, parObsName, parObsErr, parArms, time, mpOpt
##
model.defaultAE <- function() {

	default <- list()
	formPK <- myBolus_1cpt_Vk()[[1]]
	formPD <- myImmed_lin_null(formPK)[[1]]
	
	default$PFIMmodel <- c(formPK, formPD)
	default$Type <- 'AE'

	default$parModelName <- c('V','k','Alin')
	default$parModel <- c(4.5,0.5,3)
	default$parModelVar <- c(0.2,0.4,0.15)
	default$parModelVarType <- 'exp'
	default$parObsName <- c('Conc','Effect')
	default$parObsErr <- list(c(0.3,0.5), c(0,0.3))
	default$parObsTimes <- list(c(0.5,1,2,4,8,12,24),c(0.5,1,2,4,8,12))
	default$parArmsName <- c('dose','tst')
	default$parArms  <- list(c(30,5), c(10,2))
	default$ArmsName <- list('30 mg','10 mg')
	default$TimeName <- 'Time (hr)'
	default$tRange <- c(0,24)
	default$mpOpt <- list()

	# Return value: example structure
	return(default)
}

myBolus_1cpt_Vk<-function()
{
	form1<-paste("dose/V*(exp(-k*t))")
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

myImmed_lin_null<-function(fconc)
{
  form1<-paste("Alin*",fconc)
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}
