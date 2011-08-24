# Author: Martin Fink
# Date: July 2011
#
# Compare myPFIMv2 with PFIM3.2.2
#
## Input
# none
## Return list
# none
##
rm(list=ls(all=TRUE))

# Program testPFIM
myDir <- getwd()
#source(paste(myDir,"/../../tstPFIM3.2.1/PFIM3.2.r",sep=''))
source(paste(myDir,"/../../PFIM3.2.2/PFIM3.2.2/PFIM3.2.r",sep=''))
model.file="stdin.r"
file.model="model.r"
runExample <- function(n=1) {
	directory<<-paste("../../Examples/Example_",n,sep='')
	PFIM()
}

# Program myPFIMv2
names <- dir('../R','model.*.R')
for (iName in names) source(paste('../R/',iName,sep=''))
runMyExample <- function(example=1) {
	model <- model.readPFIM(paste('../../Examples/Example_',example,'/stdin.r',sep=''),example>'22')
	myfish <- model.fish(model)
	allfish <- myfish[[1]][[1]]*model$subject[1]
	if (length(myfish[[1]])>1) for (i in 2:length(myfish[[1]])) allfish <- allfish + myfish[[1]][[i]]*model$subject[i]
	sink(paste('output/out.Ex',example,'.r',sep=''))
	myout <- model.out.det_rse(model, allfish, out=T)
	sink()
	return(myout)
}

for (example in c('20', '21', '24', '25')) {
	system.time(ssEx<-runExample(example))
	system.time(stEx<-runMyExample(example))
	print(all.equal(ssEx$se, as.vector(stEx[[3]]$StdErr)))
}

example <- '21'
model <- model.readPFIM(paste('../../Examples/Example_',example,'/stdin.r',sep=''),example>'22')
a<-model.plot(model)
system.time(ssEx<-runExample(example))

# Repeat for debugging
	for (iName in names) source(paste('../R/',iName,sep=''))
	system.time(stEx<-runMyExample(example))
	print(all.equal(ssEx$se, as.vector(stEx[[3]]$StdErr)))
