#########################################################
source(paste(directory.program,dirsep,"CreateModel_PKPDdesign.r",sep=""))

create_formED(infusion_1cpt_VVmkm,turn_input_Imaxfull,dose=100,TInf=1)
# The differential equation system is created in the file model_created.r


