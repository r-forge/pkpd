#########################################################################
##					        		                                               ##
##				            INPUT FILE FOR PFIM 3.2                          ##
#########################################################################

protA<-list(c(10/60/24,1/24,6/24,0.5,1,2,3,4,5,6,7,8,9,10)) 
protB<-list(c(1,3,5,7,10)) 
protC<-list(c(1/24,6/24,0.5,1,2,3,4,5,6,7,8,9,10)) 

dosingEvents <- list(data=data.frame(var=1, time=c(1:9), val=0, method='add'))

parameters<-c("Kd","CLD","CLT","RateT","V")
#beta<-c(0.2,0.2,1,40,6)
beta<-c(0.2,0.2,1,40,6)

condinit<-expression(c(TD=0, TT=RateT/CLT*V))
# i.v. bolus: DoseD*WT/0.15, ss total target


#Name of the project
#-------------------- 

project<-"Example Escalating Dose (Within ODE) Nonlinear Monkey PK - Dosing every 2nd day"


#Name of the file containing the PK or PD model
#----------------------------------------------

file.model<-"model_created.r";


#Name of the output file for the results
#---------------------------------------

output<-"Stdout.r";

#RUN:  Evaluation (EVAL) or Optimisation (OPT) 
#-------------------------------------------------------
run<-"EVAL"

#Block diagonal Fisher information matrix (option<-1) or complete Information matrix (option<-2)
#----------------------------------------------------------
option<-1

#Number of responses
#--------------------------------------------------------------------

nr<-3

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
time.condinit<-(-2)

#Identical initial conditions in each elementary design (Yes=T, No=F)
#-------------------------------------------------------------
condinit.identical<-T

# If 'Yes', enter once the expression of the initial values of the system at the initial time
# else, enter the vectors of the initial conditions for each elementary design
# If initial values depend on parameters to be estimated, 
# enter this parameter into the expression without any quotation marks 
#---------------------------------------------------------

# Error tolerance for solving differential equations
#----------------------------------------------------

RtolEQ<-1e-8
AtolEQ<-1e-8
Hmax<-0

###### END DIFFERENTIAL EQUATION OPTION #################################


#Name of the fixed effects parameters
#-------------------------------------



#Fixed effects parameters values
#--------------------------------


#Number of occasions
#--------------------------------------------------------------------------------
n_occ<-1


#Random effect model (1) = additive  (2) = exponential 
#------------------------------------------------------------------
Trand<-2;


#Diagonal Matrix of variance for inter-subject random effects:
#---------------------------------------------------
omega<-diag(c(0.05,0.20,0.25,0.50,0.30))
#omega<-diag(c(0.0,0.0,0.0,0.0,0.0))

#Diagonal Matrix of variance for inter-occasion random effects:
#---------------------------------------------------
#gamma<-diag(c(0,0,0))

#Standard deviation of residual error (sig.inter+sig.slope*f)^2:
#------------------------------------------------------------------
sig.interA<-0
sig.slopeA<-0.2

sig.interB<-0
sig.slopeB<-0.1

sig.interC<-0
sig.slopeC<-0.1

#List of the vectors of sampling times for each elementary design 
#you can specify any sampling times for a group by writing NULL 
#ONLY if you have several responses
#-----------------------------------------------------------------
#protB<-protA
#protC<-protA 

#Vector of initial proportions or numbers of subjects for each elementary design 
#--------------------------------------------------------------
subjects<-c(6)

#Subjects input: (1) for number of subjects (2) for proportions of subjects
#---------------------------------------------------------------------------
subjects.input<-1

#If 'proportions of subjects' give the total number of samples
#-------------------------------------------------------------
#Ntot<-1000


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
parameter.associated<-list(Sex=c("V"))

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
identical.times<-F

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

subjects.opt<-T

#Vector of lower and upper admissible sampling times
#---------------------------------------------------
lowerA<-c(0) 
upperA<-c(150)

lowerB<-c(0)
upperB<-c(150)

#Minimum delay between two sampling times
#-------------------------------------------

delta.time<-0.5

#Print iteration step (Yes=T, No=F)
#---------------------------------

iter.print<-T


#Parameter for initial simplex building (%)
#------------------------------------------

simplex.parameter<-20


#Maximum iteration number
#------------------------

Max.iter<-5000


#Relative convergence tolerance
#------------------------------ 
Rctol<-1e-6



#############################
#FEDOROV-WYNN SPECIFICATION #
#############################


#Number of sampling windows
#--------------------------
nwindA<-1
nwindB<-1
nwindC<-1

#List of vector of the allowed sampling times for each sampling window
#--------------------------------------------------------------------
sampwinA<-protA
sampwinB<-protB
sampwinC<-protC

#Maximum total number of sampling times per subject
#--------------------------------------------------
nmaxptsA<-length(protA[[1]])
nmaxptsB<-length(protB[[1]])
nmaxptsC<-length(protC[[1]])

#Minimum total number of sampling times per subject
#--------------------------------------------------
nminptsA<-nmaxptsA
nminptsB<-nmaxptsB
nminptsC<-nmaxptsC

#List of vector of allowed number of points to be taken from each sampling window
#------------------------------------------------------------------------------
nsampA<-list(c(nmaxptsA))
nsampB<-list(c(nmaxptsB))
nsampC<-list(c(nmaxptsC))

############# END OF OPTIMISATION ALGORITHM OPTION ###############


############## GRAPH SPECIFICATION OPTION ###############

#graphical representation (Yes=T, No=F)
#-------------------------------------
graph.logical<-T

#Vector of Names on Y axes for each response
#---------------------------------
names.datax<-c("Time","Time","Time")

#Vector of Names on Y axes for each response
#---------------------------------
names.datay<-c("Conc free drug","Conc tot target","Conc free target")

#Controls logarithmic axes for the graphical representation.
#Values "xy", "x", or "y" produce log-log or log-x or log-y axes.
#(For standard graphic, log.logical<-F)
#--------------------------------------------------------------
log.logical<-c('y',F,F)
#log.logical<-F

#Vector of lower and upper sampling times for the graphical representation
#-------------------------------------------------------------------------
graph.infA<-c(0)
graph.supA<-c(max(protA[[1]]))
graph.infB<-c(0)
graph.supB<-c(max(protA[[1]]))
graph.infC<-c(0)
graph.supC<-c(max(protA[[1]]))

#Vector of lower and upper concentration for the graphical representation
#------------------------------------------------------------------------
y.rangeA<-c(1e-7,1e2) # default range
y.rangeB<-c(0,10000) # default range
y.rangeC<-NULL # default range
#y.range<-c(0,10)

############# END OF GRAPH SPECIFICATION OPTION ###############
