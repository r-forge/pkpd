# Author: Martin Fink
# Date: March 2011
#
# run.default <- function()
# Run default model with myPFIMv2
## Input
# none
## Return list
# none
##
rm(list=ls(all=TRUE))
names <- dir('.','model.*.R')
for (iName in names) source(iName)

runAE <- TRUE

if (runAE) {
	setup <- model.defaultAE()
	modelAE <- model.define(setup$PFIMmodel, setup$Type, setup$parModelName, setup$parModel, setup$parModelVar, setup$parModelVarType, setup$parObsName, setup$parObsErr, setup$parObsTimes, setup$parArms, setup$parArmsName, setup$ArmsName, setup$TimeName, setup$tRange, setup$parArmsIC)
	myresAE <- model.solve(modelAE)
	mysensAE <- model.sens(modelAE)
	myfishAE <- model.fish(modelAE, mysensAE)
	outAE <- model.out.det_rse(modelAE, myfishAE[[1]][[1]], out=T)

	myres2 <- model.plot(modelAE, FALSE, FALSE)
	myres2 <- model.plot(modelAE, FALSE, TRUE)
} else {
	setup <- model.defaultDE()
	modelDE <- model.define(setup$PFIMmodel, setup$Type, setup$parModelName, setup$parModel, setup$parModelVar, setup$parModelVarType, setup$parObsName, setup$parObsErr, setup$parObsTimes, setup$parArms, setup$parArmsName, setup$ArmsName, setup$TimeName, setup$tRange, setup$parArmsIC)
	myresDE <- model.solve(modelDE)
	mysensDE <- model.sens(modelDE)
	myfishDE <- model.fish(modelDE, mysensDE)
	outDE <- model.out.det_rse(modelDE, myfishDE[[1]][[1]], out=T)
}

model.plot(modelDE, T, F)
