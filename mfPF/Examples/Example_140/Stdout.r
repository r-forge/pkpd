PFIM 3.2  
 
Option: 1 

 
Project: Example Escalating Dose (Within ODE) Nonlinear Monkey PK - Dosing every 2nd day
 
Date: Mon Jan 24 00:35:50 2011
 

 
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
StartDose = 0.03 # For multi-dosing...
ndoses = 4
startt= 0                                        #days start time for dosing
tau   = 2                                        #days interdose interval of antibody

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

Error tolerance for solving differential equations system: RtolEQ = 1e-08 , AtolEQ = 1e-08 , Hmax =  0
 
Initial population design: 

 
Sample times for response: A 
                                                                               Protocol
1 c=(0.00694444444444444, 0.0416666666666667, 0.25, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  subjects                      condinit
1        6 c(TD = 0, TT = RateT/CLT * V)

 
Sample times for response: B 
            Protocol subjects                      condinit
1 c=(1, 3, 5, 7, 10)        6 c(TD = 0, TT = RateT/CLT * V)

 
Sample times for response: C 
                                                          Protocol subjects
1 c=(0.0416666666666667, 0.25, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)        6
                       condinit
1 c(TD = 0, TT = RateT/CLT * V)

 
Total number of samples: 192
 
Associated criterion value: 53.6541
 
Identical sampling times for each response: FALSE
 
Random effect model: Trand =  2
 
Variance error model response A : ( 0 + 0.2 *f)^2
Variance error model response B : ( 0 + 0.1 *f)^2
Variance error model response C : ( 0 + 0.1 *f)^2

 

Optimization step:  
 
Sampling windows for the response: A 
Window 1 : t= 0.006944444 0.04166667 0.25 0.5 1 2 3 4 5 6 7 8 9 10 
    Nb of sampling points to be taken in this window, n[ 1 ]= 14 
Maximum total number of points in one elementary protocol : 14 
Minimum total number of points in one elementary protocol : 14 

 
Sampling windows for the response: B 
Window 1 : t= 1 3 5 7 10 
    Nb of sampling points to be taken in this window, n[ 1 ]= 5 
Maximum total number of points in one elementary protocol : 5 
Minimum total number of points in one elementary protocol : 5 

 
Sampling windows for the response: C 
Window 1 : t= 0.04166667 0.25 0.5 1 2 3 4 5 6 7 8 9 10 
    Nb of sampling points to be taken in this window, n[ 1 ]= 13 
Maximum total number of points in one elementary protocol : 13 
Minimum total number of points in one elementary protocol : 13 

 

Now evaluating the Fisher Information Matrix for the 1 protocols generated 

 
**************************** OPTIMISED DESIGN *****************************
 

 
Optimised population design: 
Sample times for response: A 
                                                                                 times
1 c(0.00694444444444444, 0.0416666666666667, 0.25, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  freq Subjects                      condinit
1    1        6 c(TD = 0, TT = RateT/CLT * V)

 
Sample times for response: B 
              times freq Subjects                      condinit
1 c(1, 3, 5, 7, 10)    1        6 c(TD = 0, TT = RateT/CLT * V)

 
Sample times for response: C 
                                                            times freq Subjects
1 c(0.0416666666666667, 0.25, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)    1        6
                       condinit
1 c(TD = 0, TT = RateT/CLT * V)

 

 
Associated optimised criterion: 53.6541
 

 
******************* POPULATION FISHER INFORMATION MATRIX ******************
 
              [,1]          [,2]         [,3]         [,4]        [,5]
 [1,] 2772.1154681  95.706783614 -12.29182032 -0.146981145 0.755617141
 [2,]   95.7067890 100.689759933  -0.01042985 -0.006406578 0.392149443
 [3,]  -12.2918379  -0.010429371  19.45865238 -0.055321881 0.166443247
 [4,]   -0.1469814  -0.006406572  -0.05532188  0.006815104 0.002177600
 [5,]    0.7556177   0.392149442   0.16644320  0.002177599 0.545573181
 [6,]    0.0000000   0.000000000   0.00000000  0.000000000 0.000000000
 [7,]    0.0000000   0.000000000   0.00000000  0.000000000 0.000000000
 [8,]    0.0000000   0.000000000   0.00000000  0.000000000 0.000000000
 [9,]    0.0000000   0.000000000   0.00000000  0.000000000 0.000000000
[10,]    0.0000000   0.000000000   0.00000000  0.000000000 0.000000000
[11,]    0.0000000   0.000000000   0.00000000  0.000000000 0.000000000
[12,]    0.0000000   0.000000000   0.00000000  0.000000000 0.000000000
[13,]    0.0000000   0.000000000   0.00000000  0.000000000 0.000000000
              [,6]         [,7]         [,8]         [,9]       [,10]
 [1,] 0.000000e+00 0.000000e+00 0.000000e+00 0.0000000000  0.00000000
 [2,] 0.000000e+00 0.000000e+00 0.000000e+00 0.0000000000  0.00000000
 [3,] 0.000000e+00 0.000000e+00 0.000000e+00 0.0000000000  0.00000000
 [4,] 0.000000e+00 0.000000e+00 0.000000e+00 0.0000000000  0.00000000
 [5,] 0.000000e+00 0.000000e+00 0.000000e+00 0.0000000000  0.00000000
 [6,] 1.024617e+03 1.221305e+00 5.036302e-01 0.1152186084  0.06851492
 [7,] 1.221305e+00 1.351790e+00 3.625894e-07 0.0002189024  0.01845374
 [8,] 5.036302e-01 3.625894e-07 3.155326e+01 0.4080680730  0.08311004
 [9,] 1.152186e-01 2.189024e-04 4.080681e-01 9.9084034035  0.02276130
[10,] 6.851492e-02 1.845374e-02 8.311004e-02 0.0227613042 32.14621034
[11,] 2.354871e+01 5.775721e+00 4.906029e+00 1.2822453092  0.73701404
[12,] 5.925825e-01 5.252100e-01 6.959736e-01 0.0234252982  0.11304521
[13,] 2.757297e+01 2.145787e+01 2.122920e+01 5.0767091363  1.15651632
            [,11]        [,12]        [,13]
 [1,]    0.000000    0.0000000     0.000000
 [2,]    0.000000    0.0000000     0.000000
 [3,]    0.000000    0.0000000     0.000000
 [4,]    0.000000    0.0000000     0.000000
 [5,]    0.000000    0.0000000     0.000000
 [6,]   23.548707    0.5925825    27.572973
 [7,]    5.775721    0.5252100    21.457872
 [8,]    4.906029    0.6959736    21.229204
 [9,]    1.282245    0.0234253     5.076709
[10,]    0.737014    0.1130452     1.156516
[11,] 1934.875415   13.7012764   173.067432
[12,]   13.701276 5259.5630171   248.562660
[13,]  173.067432  248.5626597 11983.800875

 

 
************************** EXPECTED STANDARD ERRORS ************************
 
------------------------ Fixed Effects Parameters -------------------------
 
      Beta    StdError       RSE  
Kd     0.2  0.01935945  9.679723 %
CLD    0.2  0.10146426 50.732129 %
CLT    1.0  0.23015461 23.015461 %
RateT 40.0 12.27910992 30.697775 %
V      6.0  1.35903058 22.650510 %

 

 
------------------------- Variance of Inter-Subject Random Effects ----------------------
 
      Omega  StdError       RSE  
Kd     0.05 0.0312603  62.52060 %
CLD    0.20 0.8782045 439.10226 %
CLT    0.25 0.1782124  71.28494 %
RateT  0.50 0.3178147  63.56294 %
V      0.30 0.1763764  58.79213 %

 
------------------------ Standard deviation of residual error ------------------------
 
           Sigma    StdError       RSE  
sig.slopeA   0.2 0.022890939 11.445469 %
sig.slopeB   0.1 0.013795597 13.795597 %
sig.slopeC   0.1 0.009279777  9.279777 %

 
******************************* DETERMINANT ********************************
 
3.053766e+22
 
******************************** CRITERION *********************************
 
53.65408
 

 

 
Time difference of 40.453 secs
sys.self 
    0.01 

 
