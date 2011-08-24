#########################################################################
##					        		                                               ##
##				            INPUT FILE FOR PFIM 3.2                          ##
#########################################################################


#Name of the project
#-------------------- 

project<-"Example 2"


#Name of the file containing the PK or PD model
#----------------------------------------------

file.model<-"model.r";


#Name of the output file for the results
#---------------------------------------

output<-"Stdout.r";

#RUN:  Evaluation (EVAL) or Optimisation (OPT) 
#-------------------------------------------------------
run<-"EVAL"

#Block diagonal Fisher information matrix (option<-1) or complete Information matrix (option<-2)
#----------------------------------------------------------
option<-2

#Number of responses
#--------------------------------------------------------------------

nr<-1

################### MODEL OPTION ###########################

#Model form: Differential equations (DE) or analytical form (AF)
#---------------------------------------------------------------

modelform<-"DE"

###### ANALYTICAL MODEL OPTION #############################
############################################################

#Identical dose in each elementary design (Yes=T, No=F)
#-------------------------------------------------------------
#dose.identical<-T

# If 'Yes', enter the value of the dose, 
# else, enter the vector of the dose values for each elementary design
#--------------------------------------------------------------------
#dose<-c(30)

#Vector of the times intervals of each expression  
#-----------------------------------------------------------
#boundA<-list(c(0,Inf))

###### END ANALYTICAL MODEL OPTION ########################



###### DIFFERENTIAL EQUATION OPTION ##################################
######################################################################

#Initial time for which initial conditions are given
#---------------------------------------------------

time.condinit<-0

#Identical initial conditions in each elementary design (Yes=T, No=F)
#-------------------------------------------------------------
condinit.identical<-F

# If 'Yes', enter once the expression of the initial values of the system at the initial time
# else, enter the vectors of the initial conditions for each elementary design
# If initial values depend on parameters to be estimated, 
# enter this parameter into the expression without any quotation marks 
#---------------------------------------------------------

condinit<-c(expression(c(0)))


# Error tolerance for solving differential equations
#----------------------------------------------------

RtolEQ<-1e-08
AtolEQ<-1e-08
Hmax<-Inf# Default value 


###### END DIFFERENTIAL EQUATION OPTION #################################


#Name of the fixed effects parameters
#-------------------------------------

parameters<-c("ka","V","Vm","km")


#Fixed effects parameters values
#--------------------------------
beta<-c(2.72,12.2,1.0004,0.37)


#Number of occasions
#--------------------------------------------------------------------------------
n_occ<-1


#Random effect model (1) = additive  (2) = exponential 
#------------------------------------------------------------------

Trand<-2;


#Diagonal Matrix of variance for inter-subject random effects:
#---------------------------------------------------

omega<-diag(c(0.25,0.25,0.25,0.25))

#Diagonal Matrix of variance for inter-occasion random effects:
#---------------------------------------------------

gamma<-diag(c(0,0,0,0))

#Standard deviation of residual error (sig.inter+sig.slope*f)^2:
#------------------------------------------------------------------

sig.interA<-0
sig.slopeA<-0.2

#List of the vectors of sampling times for each elementary design 
#you can specify any sampling times for a group by writing NULL 
#ONLY if you have several responses
#-----------------------------------------------------------------

protA<-list(c(0.5, 2,16, 30)) 

#Vector of initial proportions or numbers of subjects for each elementary design 
#--------------------------------------------------------------

subjects<-c(30)

#Subjects input: (1) for number of subjects (2) for proportions of subjects
#---------------------------------------------------------------------------
subjects.input<-1

#If 'proportions of subjects' give the total number of samples
#-------------------------------------------------------------
Ntot<-1000





###################################################################
#                                                                 #
#                        Covariate model                          #
#                                                                 #
###################################################################

##########################################
# Covariates not changing with occasion  # 
##########################################

#Add covariate to the model  (Yes==T No==F)
#---------------------------------------------------------------------------
covariate.model<-F

#Vector of covariates
#---------------------------------------------------------------------
covariate.name<-list(c("Sex"))

#Categories for each covariate (the first category is the reference)
#-----------------------------------------------------------------------
covariate.category<-list(Sex=c("M","F"))

#Proportions of subjects in each category
#-------------------------------------------------------------------------
covariate.proportions<-list(Sex=c(0.5,0.5))

#Parameter(s) associated with each covariate
#-------------------------------------------------------------------------
parameter.associated<-list(Sex=c("ka"))

# Values of covariate parameters in covariate model 
# (values of parameters for all other categories than the reference category (for which beta=0) 
# covariate is additive on parameter if additive random effect model (Trand=1)
# covariate is additive on log parameters if exponential random effect model (Trand=2)
#-----------------------------------------------------------------------
beta.covariate<-list(Sex=list(c(log(1.2))))



#####################################
#Covariates changing with occasion  # 
#####################################


#Add covariate to the model   (Yes==T No==F)
#---------------------------------------------------------------------------
covariate_occ.model<-F

#Vector of covariates depending on the occasion
#---------------------------------------------------------------------
covariate_occ.name<-list(c("Treat"))

#Categories for each covariate (the first category is the reference)
#-----------------------------------------------------------------------
covariate_occ.category<-list(
Treat=c("A","B"))

#Sequences of values of covariates at each occasion 
#Specify as many values in each sequence as number of occasions (n_occ) for each covariate
#-------------------------------------------------------------------------------------------------------
 
covariate_occ.sequence<-list(
Treat=list(c("A","B"),c("B","A")))

#Proportions of elementary designs corresponding to each sequence of covariate values
#Specify as many values of proportion as number of sequences defined in covariate_occ.sequence for each covariate
#-----------------------------------------------------------------------------------------------------------------
covariate_occ.proportions<-list(
Treat=c(0.5,0.5))

#Parameter(s) associated with each covariate
#-------------------------------------------------------------------------
parameter_occ.associated<-list(Treat=c("Cl"))


# Values of covariate parameters in covariate model 
# (values of parameters for all other categories than the reference category (for which beta=0) 
# covariate is additive on parameter if additive random effect model (Trand=1)
# covariate is additive on log parameters if exponential random effect model (Trand=2)
#-----------------------------------------------------------------------
beta.covariate_occ<-list(Treat=list(c(log(1.1))))


#############################################
# Power and number of subjects              #
#############################################

#Type one error alpha 
#-----------------------------------------------------------------------------
alpha<-0.05

#Compute expected power for comparison test (Yes=T, No=F)
#---------------------------------------------------------------------------
compute.power<-F

#Compute the number of subjects needed for a given power for comparison test(Yes=T, No=F)
#----------------------------------------------------------------------------
compute.nni<-F

#Equivalence interval
interval_eq<-c(log(0.8),log(1.25))

#Compute expected power for equivalence test (Yes=T, No=F)
#---------------------------------------------------------------------------
compute.power_eq<-F

#Compute the number of subjects needed for a given power for equivalence test (Yes=T, No=F)
#----------------------------------------------------------------------------
compute.nni_eq<-F

#Set value the given power
#---------------------------------------------------------------------------
given.power<-0.9



############ONLY FOR OPTIMISATION ###############################

#Identical sampling times for each response
# (only if you do not have sampling times==NULL)
#-------------------------------------------------------------------------------------
#identical.times<-T

######## OPTIMISATION ALGORITHM OPTION ###############

#Character string for thoice of the optimisation algorithm: 
#	"FW" for the Fedorov-Wynn algorithm 
#	"SIMP" for the Simplex algorithm
#------------------------------------------

algo.option<-"FW"


########################
#SIMPLEX SPECIFICATION #
########################

#Optimisation of the proportions of subjects: (Yes=T, No=F)
#--------------------------------------------------------------

#subjects.opt<-T

#Vector of lower and upper admissible sampling times
#---------------------------------------------------


#lowerA<-c(0) 
#upperA<-c(24)

#lowerB<-c(0)
#upperB<-c(24)

#Minimum delay between two sampling times
#-------------------------------------------

#delta.time<-0

#Print iteration step (Yes=T, No=F)
#---------------------------------

#iter.print<-T


#Parameter for initial simplex building (%)
#------------------------------------------

#simplex.parameter<-20


#Maximum iteration number
#------------------------

#Max.iter<-5000


#Relative convergence tolerance
#------------------------------ 
#Rctol<-1e-6



#############################
#FEDOROV-WYNN SPECIFICATION #
#############################


#Number of sampling windows
#--------------------------
#nwindA<-1
#nwindB<-1


#List of vector of the allowed sampling times for each sampling window
#--------------------------------------------------------------------

#sampwinA<-list(c(0.5, 2, 5, 16, 18,30))
#sampwinB<-sampwinA


#List of vector of allowed number of points to be taken from each sampling window
#------------------------------------------------------------------------------

#nsampA<-list(c(4))
#nsampB<-list(c(4))

#Maximum total number of sampling times per subject
#--------------------------------------------------

#nmaxptsA<-4
#nmaxptsB<-4

#Minimum total number of sampling times per subject
#--------------------------------------------------

#nminptsA<-4
#nminptsB<-4
############# END OF OPTIMISATION ALGORITHM OPTION ###############






############## GRAPH SPECIFICATION OPTION ###############

#graphical representation (Yes=T, No=F)
#-------------------------------------

graph.logical<-T

#Vector of Names on Y axes for each response
#---------------------------------
names.datax<-c("Time")

#Vector of Names on Y axes for each response
#---------------------------------
names.datay<-c("Concentration")

#Controls logarithmic axes for the graphical representation.
#Values "xy", "x", or "y" produce log-log or log-x or log-y axes.
#(For standard graphic, log.logical<-F)
#--------------------------------------------------------------

#log.logical<-'y'
log.logical<-F

#Vector of lower and upper sampling times for the graphical representation
#-------------------------------------------------------------------------

graph.infA<-c(0)
graph.supA<-c(30)


#Vector of lower and upper concentration for the graphical representation
#------------------------------------------------------------------------
y.rangeA<-NULL # default range
#y.range<-c(0,10)

############# END OF GRAPH SPECIFICATION OPTION ###############













