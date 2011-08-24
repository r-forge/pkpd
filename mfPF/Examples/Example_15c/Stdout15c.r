PFIM 3.2  
 
Option: 1 

 
Project: Example Escalating Dose (Within ODE) Nonlinear Monkey PK - Dosing every 2nd day - only PK samples
 
Date: Thu Jan 13 12:11:54 2011
 

 
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
  subjects               condinit
1      0.5 c(0, RateT/CLT * V, 0)
2      0.5 c(0, RateT/CLT * V, 1)

 
Sample times for response: B 
  Protocol subjects               condinit
1    c=(2)      0.5 c(0, RateT/CLT * V, 0)
2    c=(2)      0.5 c(0, RateT/CLT * V, 1)

 
Total number of samples: 104
 
Associated criterion value: 41.3705
 
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
                                                                   times
1 c(0.0416666666666667, 1, 2, 2.04166666666667, 3, 4, 5, 6, 7, 8, 9, 10)
2 c(0.0416666666666667, 1, 2, 2.04166666666667, 3, 4, 5, 6, 7, 8, 9, 10)
       freq Subjects               condinit
1 0.2066086 1.652869 c(0, RateT/CLT * V, 0)
2 0.7933914 6.347131 c(0, RateT/CLT * V, 1)

 
Sample times for response: B 
  times      freq Subjects               condinit
1     2 0.2066086 1.652869 c(0, RateT/CLT * V, 0)
2     2 0.7933914 6.347131 c(0, RateT/CLT * V, 1)

 

 
Associated optimised criterion: 41.3705
 

 
******************* POPULATION FISHER INFORMATION MATRIX ******************
 
             [,1]         [,2]       [,3]         [,4]         [,5]
 [1,] 10.50367685   7.56911247 -2.2470997 -0.040027068 -0.114592218
 [2,]  7.56911247 164.41091330  9.7356654  0.085112529  3.012228531
 [3,] -2.24709965   9.73566545 19.7051289 -0.143928891  0.257498004
 [4,] -0.04002707   0.08511253 -0.1439289  0.008173337  0.005239479
 [5,] -0.11459222   3.01222853  0.2574980  0.005239479  2.095422065
 [6,]  0.00000000   0.00000000  0.0000000  0.000000000  0.000000000
 [7,]  0.00000000   0.00000000  0.0000000  0.000000000  0.000000000
 [8,]  0.00000000   0.00000000  0.0000000  0.000000000  0.000000000
 [9,]  0.00000000   0.00000000  0.0000000  0.000000000  0.000000000
[10,]  0.00000000   0.00000000  0.0000000  0.000000000  0.000000000
[11,]  0.00000000   0.00000000  0.0000000  0.000000000  0.000000000
[12,]  0.00000000   0.00000000  0.0000000  0.000000000  0.000000000
[13,]  0.00000000   0.00000000  0.0000000  0.000000000  0.000000000
              [,6]        [,7]       [,8]        [,9]       [,10]        [,11]
 [1,]    0.0000000  0.00000000  0.0000000  0.00000000   0.0000000 0.000000e+00
 [2,]    0.0000000  0.00000000  0.0000000  0.00000000   0.0000000 0.000000e+00
 [3,]    0.0000000  0.00000000  0.0000000  0.00000000   0.0000000 0.000000e+00
 [4,]    0.0000000  0.00000000  0.0000000  0.00000000   0.0000000 0.000000e+00
 [5,]    0.0000000  0.00000000  0.0000000  0.00000000   0.0000000 0.000000e+00
 [6,]  139.0577084  0.72478454  1.5902075  0.80742150   0.1487752 1.122120e+03
 [7,]    0.7247845  5.49252351  0.3343465  0.07242571   0.9039779 1.905176e+00
 [8,]    1.5902075  0.33434653 25.3607282  2.36548489   0.1497185 2.845157e+01
 [9,]    0.8074215  0.07242571  2.3654849 10.77295107   0.1024275 7.431555e+00
[10,]    0.1487752  0.90397788  0.1497185  0.10242755 355.6649751 4.024929e+00
[11,] 1122.1199185  1.90517561 28.4515683  7.43155517   4.0249293 1.233799e+05
[12,]   56.8651195 17.31250250 17.7847958  4.46078501  14.2785519 4.364425e+03
[13,]   35.9251935  4.05001864  4.9510201  0.24974468   8.5248423 3.434894e+01
            [,12]       [,13]
 [1,]    0.000000   0.0000000
 [2,]    0.000000   0.0000000
 [3,]    0.000000   0.0000000
 [4,]    0.000000   0.0000000
 [5,]    0.000000   0.0000000
 [6,]   56.865120  35.9251935
 [7,]   17.312502   4.0500186
 [8,]   17.784796   4.9510201
 [9,]    4.460785   0.2497447
[10,]   14.278552   8.5248423
[11,] 4364.425193  34.3489416
[12,] 3009.418365  89.3489371
[13,]   89.348937 125.3215303

 

 
************************** EXPECTED STANDARD ERRORS ************************
 
------------------------ Fixed Effects Parameters -------------------------
 
      Beta   StdError      RSE  
Kd     2.0  0.3305837 16.52919 %
CLD    0.2  0.0837738 41.88690 %
CLT    1.0  0.2568230 25.68230 %
RateT 40.0 12.3709347 30.92734 %
V      6.0  0.7013781 11.68963 %

 

 
------------------------- Variance of Inter-Subject Random Effects ----------------------
 
      Omega   StdError       RSE  
Kd     0.05 0.09181679 183.63357 %
CLD    0.20 0.43500699 217.50350 %
CLT    0.25 0.20167158  80.66863 %
RateT  0.50 0.30795085  61.59017 %
V      0.10 0.05307818  53.07818 %

 
------------------------ Standard deviation of residual error ------------------------
 
           Sigma    StdError       RSE  
sig.interA  0.01 0.003041859 30.418587 %
sig.slopeA  0.20 0.019080828  9.540414 %
sig.slopeB  0.10 0.095283619 95.283619 %

 
******************************* DETERMINANT ********************************
 
1.039870e+21
 
******************************** CRITERION *********************************
 
41.37049
 

 

 
Time difference of 14.375 secs
sys.self 
    0.17 

 
