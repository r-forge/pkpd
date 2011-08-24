PFIM 3.2  
 
Option: 1 

 
Project: Example Escalating Dose (Within ODE) Nonlinear Monkey PK - Dosing every 2nd day - only PK samples
 
Date: Mon Jan 24 00:46:56 2011
 

 
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
ndoses = 4
startt= 0                                        #days start time for dosing

TD <- y[1]
TT <- y[2]
DoseSwitch <- y[3]
#TT <- RateT/CLT*V

#DoseD = if time<0 then 0                #mg/kg dose of antibody
#        else StartDose*(int(time)+1) #increment dose each day
#input = if time>ndoses*tau-DT then 0        #input dose in nmoles
#        else pulse(DoseD*WT/0.15,startt+DT,tau)
#TDose'=input
#INIT TDose=0
BolusDuration=0.01 #  14.4/60/24
input = 0
Day <- floor(t) + 1
if (DoseSwitch==0) {
	tau   = 100        #days interdose interval of antibody
	StartDose = 10 # For multi-dosing...
} else {
	tau   = 2        #days interdose interval of antibody
	StartDose = 0.06 # For multi-dosing...
}

if ((t-startt)%%tau < BolusDuration & t>0 & t<(ndoses*tau+0.5)) input <- StartDose*Day^2/BolusDuration*WT/0.15

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
return(list(c(dTD,dTT,0),c(CFD,CTT)))
}

 
Initial Conditions at time -2 : 
 
0 RateT/CLT * V 0 

Error tolerance for solving differential equations system: RtolEQ = 1e-08 , AtolEQ = 1e-08 , Hmax =  0.01
 
Initial population design: 

 
Sample times for response: A 
                                                                 Protocol
1 c=(0.0416666666666667, 1, 2, 2.04166666666667, 3, 4, 5, 6, 7, 8, 9, 10)
2 c=(0.0416666666666667, 1, 2, 2.04166666666667, 3, 4, 5, 6, 7, 8, 9, 10)
  subjects                              condinit
1        8 c(TD = 0, TT = RateT/CLT * V, TP = 0)
2        8 c(TD = 0, TT = RateT/CLT * V, TP = 1)

 
Sample times for response: B 
  Protocol subjects                              condinit
1    c=(2)        8 c(TD = 0, TT = RateT/CLT * V, TP = 0)
2    c=(2)        8 c(TD = 0, TT = RateT/CLT * V, TP = 1)

 
Total number of samples: 104
 
Associated criterion value: 0.0631
 
Identical sampling times for each response: FALSE
 
Random effect model: Trand =  2
 
Variance error model response A : ( 0.01 + 0.2 *f)^2
Variance error model response B : ( 0 + 0.1 *f)^2

 

Optimization step:  
 
Sampling windows for the response: A 
Window 1 : t= 0.04166667 1 2 2.041667 3 4 5 6 7 8 9 10 
    Nb of sampling points to be taken in this window, n[ 1 ]= 12 
Maximum total number of points in one elementary protocol : 12 
Minimum total number of points in one elementary protocol : 12 

 
Sampling windows for the response: B 
Window 1 : t= 2 
    Nb of sampling points to be taken in this window, n[ 1 ]= 1 
Maximum total number of points in one elementary protocol : 1 
Minimum total number of points in one elementary protocol : 1 

 

Now evaluating the Fisher Information Matrix for the 2 protocols generated 

 
**************************** OPTIMISED DESIGN *****************************
 

 
Optimised population design: 
Sample times for response: A 
                                                                   times freq
1 c(0.0416666666666667, 1, 2, 2.04166666666667, 3, 4, 5, 6, 7, 8, 9, 10)    1
  Subjects                              condinit
1        8 c(TD = 0, TT = RateT/CLT * V, TP = 1)

 
Sample times for response: B 
  times freq Subjects                              condinit
1     2    1        8 c(TD = 0, TT = RateT/CLT * V, TP = 1)

 

 
Associated optimised criterion: 40.4132
 

 
******************* POPULATION FISHER INFORMATION MATRIX ******************
 
              [,1]         [,2]          [,3]          [,4]          [,5]
 [1,] 301.54891148 7.089095e-02 -56.483706069 -7.427386e-01  0.1313937366
 [2,]   0.07089094 9.999999e+02   0.007855678  7.163852e-05  0.0005477782
 [3,] -56.48370328 7.855677e-03  25.882352740 -6.007483e-02 -0.2602479696
 [4,]  -0.74273857 7.163850e-05  -0.060074831  9.279355e-03  0.0018132692
 [5,]   0.13139362 5.477785e-04  -0.260247970  1.813269e-03  2.0477356926
 [6,]   0.00000000 0.000000e+00   0.000000000  0.000000e+00  0.0000000000
 [7,]   0.00000000 0.000000e+00   0.000000000  0.000000e+00  0.0000000000
 [8,]   0.00000000 0.000000e+00   0.000000000  0.000000e+00  0.0000000000
 [9,]   0.00000000 0.000000e+00   0.000000000  0.000000e+00  0.0000000000
[10,]   0.00000000 0.000000e+00   0.000000000  0.000000e+00  0.0000000000
[11,]   0.00000000 0.000000e+00   0.000000000  0.000000e+00  0.0000000000
[12,]   0.00000000 0.000000e+00   0.000000000  0.000000e+00  0.0000000000
[13,]   0.00000000 0.000000e+00   0.000000000  0.000000e+00  0.0000000000
              [,6]         [,7]         [,8]         [,9]        [,10]
 [1,] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
 [2,] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
 [3,] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
 [4,] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
 [5,] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
 [6,] 9.093175e+00 6.076202e-07 7.976022e+00 2.206642e+00 1.553787e-03
 [7,] 6.076202e-07 9.999998e+01 4.734781e-07 5.928277e-08 9.971969e-08
 [8,] 7.976022e+00 4.734781e-07 4.186851e+01 3.608985e-01 1.523903e-01
 [9,] 2.206642e+00 5.928277e-08 3.608985e-01 1.377703e+01 1.183660e-02
[10,] 1.553787e-03 9.971969e-08 1.523903e-01 1.183660e-02 3.396509e+02
[11,] 1.019823e+02 2.675115e-06 4.134322e+01 1.089892e+01 7.778502e+00
[12,] 1.834293e+01 4.485847e-06 1.000791e+01 2.750642e+00 5.010341e+00
[13,] 2.692647e+00 1.411553e-06 1.344105e+01 7.340251e-02 4.620257e+01
             [,11]        [,12]        [,13]
 [1,] 0.000000e+00 0.000000e+00 0.000000e+00
 [2,] 0.000000e+00 0.000000e+00 0.000000e+00
 [3,] 0.000000e+00 0.000000e+00 0.000000e+00
 [4,] 0.000000e+00 0.000000e+00 0.000000e+00
 [5,] 0.000000e+00 0.000000e+00 0.000000e+00
 [6,] 1.019823e+02 1.834293e+01 2.692647e+00
 [7,] 2.675115e-06 4.485847e-06 1.411553e-06
 [8,] 4.134322e+01 1.000791e+01 1.344105e+01
 [9,] 1.089892e+01 2.750642e+00 7.340251e-02
[10,] 7.778502e+00 5.010341e+00 4.620257e+01
[11,] 5.186123e+03 1.194215e+03 1.420107e+01
[12,] 1.194215e+03 2.256184e+03 6.300406e+00
[13,] 1.420107e+01 6.300406e+00 2.569229e+01

 

 
************************** EXPECTED STANDARD ERRORS ************************
 
------------------------ Fixed Effects Parameters -------------------------
 
      Beta    StdError      RSE  
Kd     0.2  0.10276713 51.38357 %
CLD    0.2  0.03162278 15.81139 %
CLT    1.0  0.31685447 31.68545 %
RateT 40.0 14.34966161 35.87415 %
V      6.0  0.69947660 11.65794 %

 

 
------------------------- Variance of Inter-Subject Random Effects ----------------------
 
      Omega   StdError       RSE  
Kd     0.05 0.42402155 848.04311 %
CLD    0.20 0.10000001  50.00001 %
CLT    0.25 0.18967528  75.87011 %
RateT  0.50 0.27625424  55.25085 %
V      0.10 0.06455522  64.55522 %

 
------------------------ Standard deviation of residual error ------------------------
 
           Sigma   StdError       RSE  
sig.interA  0.01 0.01684335 168.43353 %
sig.slopeA  0.20 0.02249745  11.24872 %
sig.slopeB  0.10 0.25735543 257.35543 %

 
******************************* DETERMINANT ********************************
 
7.670177e+20
 
******************************** CRITERION *********************************
 
40.41322
 

 

 
Time difference of 1.5276 mins
sys.self 
    0.06 

 
