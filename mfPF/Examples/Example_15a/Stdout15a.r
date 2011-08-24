PFIM 3.2  
 
Option: 1 

 
Project: Example Escalating Dose (Within ODE) Nonlinear Monkey PK - Dosing every 2nd day - only PK samples
 
Date: Thu Jan 13 11:51:54 2011
 

 
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
1 0.3333333 c(0, RateT/CLT * V, 0)
2 0.6666667 c(0, RateT/CLT * V, 1)

 
Sample times for response: B 
  Protocol  subjects               condinit
1    c=(2) 0.3333333 c(0, RateT/CLT * V, 0)
2    c=(2) 0.6666667 c(0, RateT/CLT * V, 1)

 
Total number of samples: 78
 
Associated criterion value: 38.4546
 
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
1 0.7992279 4.795368 c(0, RateT/CLT * V, 1)
2 0.2007721 1.204632 c(0, RateT/CLT * V, 0)

 
Sample times for response: B 
  times      freq Subjects               condinit
1     2 0.7992279 4.795368 c(0, RateT/CLT * V, 1)
2     2 0.2007721 1.204632 c(0, RateT/CLT * V, 0)

 

 
Associated optimised criterion: 38.4546
 

 
******************* POPULATION FISHER INFORMATION MATRIX ******************
 
             [,1]         [,2]          [,3]         [,4]        [,5]
 [1,] 217.7291775  21.47731612 -31.014353573 -0.435521178 0.894005679
 [2,]  21.4773161 125.86941417   7.963740680  0.060421289 2.008292508
 [3,] -31.0143536   7.96374068  17.621050867 -0.072985437 0.001133466
 [4,]  -0.4355212   0.06042129  -0.072985437  0.006579568 0.001743494
 [5,]   0.8940057   2.00829251   0.001133466  0.001743494 1.585587330
 [6,]   0.0000000   0.00000000   0.000000000  0.000000000 0.000000000
 [7,]   0.0000000   0.00000000   0.000000000  0.000000000 0.000000000
 [8,]   0.0000000   0.00000000   0.000000000  0.000000000 0.000000000
 [9,]   0.0000000   0.00000000   0.000000000  0.000000000 0.000000000
[10,]   0.0000000   0.00000000   0.000000000  0.000000000 0.000000000
[11,]   0.0000000   0.00000000   0.000000000  0.000000000 0.000000000
[12,]   0.0000000   0.00000000   0.000000000  0.000000000 0.000000000
[13,]   0.0000000   0.00000000   0.000000000  0.000000000 0.000000000
              [,6]        [,7]        [,8]       [,9]        [,10]        [,11]
 [1,]   0.00000000  0.00000000  0.00000000 0.00000000   0.00000000 0.000000e+00
 [2,]   0.00000000  0.00000000  0.00000000 0.00000000   0.00000000 0.000000e+00
 [3,]   0.00000000  0.00000000  0.00000000 0.00000000   0.00000000 0.000000e+00
 [4,]   0.00000000  0.00000000  0.00000000 0.00000000   0.00000000 0.000000e+00
 [5,]   0.00000000  0.00000000  0.00000000 0.00000000   0.00000000 0.000000e+00
 [6,]   7.90860943  0.07770295  4.01058416 1.26539548   0.12221203 1.498355e+02
 [7,]   0.07770295  4.06056222  0.26805438 0.05084532   0.56491659 7.433188e-01
 [8,]   4.01058416  0.26805438 27.57402708 1.15262154   0.01575584 2.596914e+01
 [9,]   1.26539548  0.05084532  1.15262154 9.35927185   0.01479160 6.762288e+00
[10,]   0.12221203  0.56491659  0.01575584 0.01479160 271.54246086 3.191665e+00
[11,] 149.83551710  0.74331880 25.96914103 6.76228789   3.19166529 2.544211e+05
[12,]  11.65095206 12.60183947  7.98718600 1.77308118   8.22053118 1.563283e+03
[13,]   7.63549433  5.77124467  3.12855838 0.19650835   8.35771460 2.689971e+01
            [,12]       [,13]
 [1,]    0.000000   0.0000000
 [2,]    0.000000   0.0000000
 [3,]    0.000000   0.0000000
 [4,]    0.000000   0.0000000
 [5,]    0.000000   0.0000000
 [6,]   11.650952   7.6354943
 [7,]   12.601839   5.7712447
 [8,]    7.987186   3.1285584
 [9,]    1.773081   0.1965083
[10,]    8.220531   8.3577146
[11,] 1563.282617  26.8997096
[12,] 2030.196594  87.5840524
[13,]   87.584052 217.0419477

 

 
************************** EXPECTED STANDARD ERRORS ************************
 
------------------------ Fixed Effects Parameters -------------------------
 
      Beta    StdError      RSE  
Kd     0.2  0.10181527 50.90763 %
CLD    0.2  0.09892295 49.46148 %
CLT    1.0  0.34257254 34.25725 %
RateT 40.0 15.97486148 39.93715 %
V      6.0  0.80284188 13.38070 %

 

 
------------------------- Variance of Inter-Subject Random Effects ----------------------
 
      Omega  StdError       RSE  
Kd     0.05 0.3822193 764.43855 %
CLD    0.20 0.5098191 254.90956 %
CLT    0.25 0.1980983  79.23931 %
RateT  0.50 0.3307777  66.15554 %
V      0.10 0.0607257  60.72570 %

 
------------------------ Standard deviation of residual error ------------------------
 
           Sigma    StdError      RSE  
sig.interA  0.01 0.001998373 19.98373 %
sig.slopeA  0.20 0.022647306 11.32365 %
sig.slopeB  0.10 0.070806270 70.80627 %

 
******************************* DETERMINANT ********************************
 
4.020882e+20
 
******************************** CRITERION *********************************
 
38.45456
 

 

 
Time difference of 14.125 secs
sys.self 
    0.08 

 
