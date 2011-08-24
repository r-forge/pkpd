# Author: Martin Fink
# Date: July 2011
#
# Run diverse models and compare the sensitivity arrays
#
## Input
# none
## Return list
# none
##
rm(list=ls(all=TRUE))

# Program myPFIMv2
names <- dir('../R','model.*.R')
for (iName in names) source(paste('../R/',iName,sep=''))

runMyModel <- function(modelJ) {
	# Derive sensitivities with 3 different methods
	print(system.time(mHessJ <- model.sensNum(modelJ, useHess=T)))
	print(system.time(mJac <- model.sensNum(modelJ)))
	print(system.time(mPDJ <- model.sens(modelJ)))
	print(all.equal(mPDJ,mHessJ))
	print(all.equal(mPDJ,mJac))
	return()
}

for (example in c('20','21','24','25')) {
	cat('\nExample',example,'\n')
	modelJ <- model.readPFIM(paste('../../Examples/Example_',example,'/stdin.r',sep=''),example>'22')
	runMyModel(modelJ)
}

for (idx in 1:9) {
	cat('\nDefault model',idx,'\n')
	modelJ <- model.define(switch(idx, model.defaultAE(), model.defaultAEfun(), model.defaultAEfun2(), model.defaultAEfun3(), model.defaultDE(), model.defaultDE2(), model.defaultDEjac(), model.defaultDEjac2(), model.defaultDEjac3()))
	model.plot(modelJ)
	runMyModel(modelJ)
}

for (iName in names) source(paste('../R/',iName,sep=''))
idx=3
modelJ <- model.define(switch(idx, model.defaultAE(), model.defaultAEfun(), model.defaultAEfun2(), model.defaultAEfun3(), model.defaultDE(), model.defaultDE2(), model.defaultDEjac(), model.defaultDEjac2(), model.defaultDEjac3()))
#res <- model.sens(modelJ, modelJ$fullTime)
#a <- model.plot(modelJ, res=res)
a <- model.plot(modelJ)
a <- model.plot(modelJ,F,T)
