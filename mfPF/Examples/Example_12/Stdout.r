PFIM 3.2  
 
Option: 1 

 
Project: Example Escalating Dose Nonlinear Monkey PK
 
Date: Wed Jan 12 13:51:08 2011
 

 
**************************** INPUT SUMMARY ********************************
 
Differential Equations form of the model:  
 
function(t,y,p){

Kd  <- p[1]
CLD <- p[2]
CLT <- p[3]
RateT <- p[4]
V <- p[5]
#Kd    = 0.2                #nM equilibrium dissociation constant
#CLD   = 0.2                #(L/d) elimination rate constant for the drug
#CLT   = 1                #(L/d) elimination rate constant for the binding target
CLC   = CLT                #(L/d) elimination rate constant for the TD complex
#RateT = 40				#(nmoles/d) input rate, or expression, of the target
#V     = 6                #(L) volume in which the binding reaction takes place

WT    = 70                                        #kg
StartDose = 5#0.03 # For multi-dosing...
ndoses = 10
startt= 0                                        #days start time for dosing
tau   = 100                                      #days interdose interval of antibody

#DoseD = if time<0 then 0                #mg/kg dose of antibody
#        else StartDose*(int(time)+1) #increment dose each day
#input = if time>ndoses*tau-DT then 0        #input dose in nmoles
#        else pulse(DoseD*WT/0.15,startt+DT,tau)
#TDose'=input
#INIT TDose=0
BolusDuration=0.01 #  14.4/60/24
input = 0
Day <- floor(t) + 1
if ((t-startt)%%tau < BolusDuration & t>0 & t<(ndoses*tau+0.5)) input <- StartDose*Day^2/BolusDuration*WT/0.15

TD <- y[1]
TT <- y[2]
#TT <- RateT/CLT*V

C =((Kd*V+TD+TT)-sqrt((Kd*V+TD+TT)^2 -4*TD*TT))/2
D = TD-C        #free drug: total drug minus complex
T = TT-C        #free target: total target minus complex

dTD <- input - D*CLD/V - C*CLC/V	# Drug
dTT <- RateT - T*CLT/V - C*CLC/V	# Drug-Target Complex

#where
#INIT TD = 0 - or i.v. bolus: DoseD*WT/0.15
#INIT TT = RateT/CLT*V        #initial target is rate in divided by elimination

CFD = 0.15*D/V        #µg/mL total antibody is free plus complex
CTT = 190 *TT/V    #ng/mL free plus complex for total target
CFT = 190 *T/V        #ng/mL free target

# if(t<=1){
#	dpk1<-(100/(1*V))+(-Vm)*pk[1]/(km*V+pk[1])}
# else{
#	dpk1<-(-Vm)*pk[1]/(km*V+pk[1])}
return(list(c(dTD,dTT),c(CFD,CTT,CFT)))
}

 
Initial Conditions at time -2 : 
 
0 RateT/CLT * V 

Error tolerance for solving differential equations system: RtolEQ = 1e-08 , AtolEQ = 1e-08 , Hmax =  0.01
 
Initial population design: 

 
Sample times for response: A 
                                                                                     Protocol
1 c=(0.00694444444444444, 0.0416666666666667, 0.25, 0.5, 2, 4, 7, 14, 21, 35, 49, 63, 77, 91)
  subjects            condinit
1        1 c(0, RateT/CLT * V)

 
Sample times for response: B 
                                Protocol subjects            condinit
1 c=(0.00694444444444444, 7, 35, 77, 91)        1 c(0, RateT/CLT * V)

 
Sample times for response: C 
                                      Protocol subjects            condinit
1 c=(0.5, 2, 4, 7, 14, 21, 35, 49, 63, 77, 91)        1 c(0, RateT/CLT * V)

 
Total number of samples: 180
 
Associated criterion value: 86.5002
 
Identical sampling times for each response: FALSE
 
Random effect model: Trand =  2
 
Variance error model response A : ( 0 + 0.2 *f)^2
Variance error model response B : ( 0 + 0.1 *f)^2
Variance error model response C : ( 0 + 0.1 *f)^2

 

Optimization step:  
 
Sampling windows for the response: A 
Window 1 : t= 0.006944444 0.04166667 0.25 0.5 2 4 7 14 21 35 49 63 77 91 
    Nb of sampling points to be taken in this window, n[ 1 ]= 14 
Maximum total number of points in one elementary protocol : 14 
Minimum total number of points in one elementary protocol : 14 

 
Sampling windows for the response: B 
Window 1 : t= 0.006944444 7 35 77 91 
    Nb of sampling points to be taken in this window, n[ 1 ]= 5 
Maximum total number of points in one elementary protocol : 5 
Minimum total number of points in one elementary protocol : 5 

 
Sampling windows for the response: C 
Window 1 : t= 0.5 2 4 7 14 21 35 49 63 77 91 
    Nb of sampling points to be taken in this window, n[ 1 ]= 11 
Maximum total number of points in one elementary protocol : 11 
Minimum total number of points in one elementary protocol : 11 

 

Now evaluating the Fisher Information Matrix for the 1 protocols generated 

 
**************************** OPTIMISED DESIGN *****************************
 

 
Optimised population design: 
Sample times for response: A 
                                                                                       times
1 c(0.00694444444444444, 0.0416666666666667, 0.25, 0.5, 2, 4, 7, 14, 21, 35, 49, 63, 77, 91)
  freq Subjects            condinit
1    1        6 c(0, RateT/CLT * V)

 
Sample times for response: B 
                                  times freq Subjects            condinit
1 c(0.00694444444444444, 7, 35, 77, 91)    1        6 c(0, RateT/CLT * V)

 
Sample times for response: C 
                                        times freq Subjects            condinit
1 c(0.5, 2, 4, 7, 14, 21, 35, 49, 63, 77, 91)    1        6 c(0, RateT/CLT * V)

 

 
Associated optimised criterion: 86.5002
 

 
******************* POPULATION FISHER INFORMATION MATRIX ******************
 
               [,1]        [,2]        [,3]          [,4]          [,5]
 [1,] 2629.81159414  55.7283673 -4.05427549 -0.0172424339  0.9830559553
 [2,]   55.72836669 653.4425081  8.43340244  0.1052126578  0.1154422028
 [3,]   -4.05427539   8.4334024 22.99441745 -0.0109876754 -0.0356955788
 [4,]   -0.01724243   0.1052127 -0.01098768  0.0073658756 -0.0004280675
 [5,]    0.98305596   0.1154422 -0.03569558 -0.0004280675  0.5487736763
 [6,]    0.00000000   0.0000000  0.00000000  0.0000000000  0.0000000000
 [7,]    0.00000000   0.0000000  0.00000000  0.0000000000  0.0000000000
 [8,]    0.00000000   0.0000000  0.00000000  0.0000000000  0.0000000000
 [9,]    0.00000000   0.0000000  0.00000000  0.0000000000  0.0000000000
[10,]    0.00000000   0.0000000  0.00000000  0.0000000000  0.0000000000
[11,]    0.00000000   0.0000000  0.00000000  0.0000000000  0.0000000000
[12,]    0.00000000   0.0000000  0.00000000  0.0000000000  0.0000000000
[13,]    0.00000000   0.0000000  0.00000000  0.0000000000  0.0000000000
              [,6]         [,7]         [,8]         [,9]        [,10]
 [1,] 0.000000e+00 0.000000e+00  0.000000000 0.000000e+00 0.000000e+00
 [2,] 0.000000e+00 0.000000e+00  0.000000000 0.000000e+00 0.000000e+00
 [3,] 0.000000e+00 0.000000e+00  0.000000000 0.000000e+00 0.000000e+00
 [4,] 0.000000e+00 0.000000e+00  0.000000000 0.000000e+00 0.000000e+00
 [5,] 0.000000e+00 0.000000e+00  0.000000000 0.000000e+00 0.000000e+00
 [6,] 9.221212e+02 4.140868e-01  0.054790498 1.585608e-03 1.159679e-01
 [7,] 4.140868e-01 5.693161e+01  0.237074255 5.903842e-02 1.599228e-03
 [8,] 5.479050e-02 2.370743e-01 44.061936166 1.609720e-02 3.822523e-03
 [9,] 1.585608e-03 5.903842e-02  0.016097201 1.157464e+01 8.795607e-04
[10,] 1.159679e-01 1.599228e-03  0.003822523 8.795607e-04 3.252448e+01
[11,] 3.616876e+01 4.447048e+00  1.243392508 2.878393e-01 9.923008e-01
[12,] 1.490622e+01 5.596377e-04  0.834353640 5.535119e-03 3.943144e-03
[13,] 3.991773e+01 2.255645e+01  5.126360345 1.202897e+00 2.728548e-01
             [,11]        [,12]        [,13]
 [1,]    0.0000000 0.000000e+00    0.0000000
 [2,]    0.0000000 0.000000e+00    0.0000000
 [3,]    0.0000000 0.000000e+00    0.0000000
 [4,]    0.0000000 0.000000e+00    0.0000000
 [5,]    0.0000000 0.000000e+00    0.0000000
 [6,]   36.1687593 1.490622e+01   39.9177301
 [7,]    4.4470482 5.596377e-04   22.5564501
 [8,]    1.2433925 8.343536e-01    5.1263603
 [9,]    0.2878393 5.535119e-03    1.2028968
[10,]    0.9923008 3.943144e-03    0.2728548
[11,] 2358.1041298 3.556348e+00  224.3410691
[12,]    3.5563484 5.041032e+03  270.5049138
[13,]  224.3410691 2.705049e+02 9336.2946224

 

 
************************** EXPECTED STANDARD ERRORS ************************
 
------------------------ Fixed Effects Parameters -------------------------
 
      Beta    StdError       RSE  
Kd     0.2  0.01952820  9.764102 %
CLD    0.2  0.03930048 19.650242 %
CLT    1.0  0.20918117 20.918117 %
RateT 40.0 11.67092382 29.177310 %
V      6.0  1.35047699 22.507950 %

 

 
------------------------- Variance of Inter-Subject Random Effects ----------------------
 
      Omega   StdError      RSE  
Kd     0.05 0.03294425 65.88850 %
CLD    0.20 0.13260561 66.30280 %
CLT    0.25 0.15065709 60.26284 %
RateT  0.50 0.29393468 58.78694 %
V      0.30 0.17534676 58.44892 %

 
------------------------ Standard deviation of residual error ------------------------
 
           Sigma   StdError      RSE  
sig.slopeA   0.2 0.02062387 10.31194 %
sig.slopeB   0.1 0.01409574 14.09574 %
sig.slopeC   0.1 0.01037519 10.37519 %

 
******************************* DETERMINANT ********************************
 
1.517820e+25
 
******************************** CRITERION *********************************
 
86.50017
 

 

 
Time difference of 53.203 secs
sys.self 
    0.11 

 
