#########################################################
source(paste(directory.program,dirsep,"CreateModel_PKPDdesign.r",sep=""))

create_formED(infusion_1cpt_VVmkm,immed_lin_null,dose=100,TInf=1)
#The differential equation system is created in the file model_created.r


