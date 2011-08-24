PFIM 3.2  
 
Option: 1 

 
Project: Example Escalating Dose (Within ODE) Nonlinear Monkey PK - Dosing every 2nd day - only PK samples
 
Date: Mon Jan 24 20:20:26 2011
 

 
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
 
Associated criterion value: 61.7404
 
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
1 0.7946225  6.35698 c(0, RateT/CLT * V, 1)
2 0.2053775  1.64302 c(0, RateT/CLT * V, 0)

 
Sample times for response: B 
  times      freq Subjects               condinit
1     2 0.7946225  6.35698 c(0, RateT/CLT * V, 1)
2     2 0.2053775  1.64302 c(0, RateT/CLT * V, 0)

 

 
Associated optimised criterion: 61.7404
 

 
******************* POPULATION FISHER INFORMATION MATRIX ******************
 
              [,1]        [,2]          [,3]         [,4]         [,5]
 [1,] 11124.200062 114.3623283 -241.03877222 -3.353203633 11.729464658
 [2,]   114.362328 165.6753788   11.36353287  0.090590306  2.724733197
 [3,]  -241.038772  11.3635329   23.79923513 -0.092798480 -0.052923984
 [4,]    -3.353204   0.0905903   -0.09279848  0.008845796  0.001570611
 [5,]    11.729465   2.7247332   -0.05292398  0.001570611  2.117795929
 [6,]     0.000000   0.0000000    0.00000000  0.000000000  0.000000000
 [7,]     0.000000   0.0000000    0.00000000  0.000000000  0.000000000
 [8,]     0.000000   0.0000000    0.00000000  0.000000000  0.000000000
 [9,]     0.000000   0.0000000    0.00000000  0.000000000  0.000000000
[10,]     0.000000   0.0000000    0.00000000  0.000000000  0.000000000
[11,]     0.000000   0.0000000    0.00000000  0.000000000  0.000000000
[12,]     0.000000   0.0000000    0.00000000  0.000000000  0.000000000
[13,]     0.000000   0.0000000    0.00000000  0.000000000  0.000000000
              [,6]        [,7]        [,8]        [,9]        [,10]
 [1,]   0.00000000  0.00000000  0.00000000  0.00000000   0.00000000
 [2,]   0.00000000  0.00000000  0.00000000  0.00000000   0.00000000
 [3,]   0.00000000  0.00000000  0.00000000  0.00000000   0.00000000
 [4,]   0.00000000  0.00000000  0.00000000  0.00000000   0.00000000
 [5,]   0.00000000  0.00000000  0.00000000  0.00000000   0.00000000
 [6,]   1.55730606  0.01682305  1.82715683  0.56577883   0.15865358
 [7,]   0.01682305  5.48765550  0.38219777  0.06865550   0.78370189
 [8,]   1.82715683  0.38219777 37.87300592  1.51028222   0.02699149
 [9,]   0.56577883  0.06865550  1.51028222 12.70380754   0.00918335
[10,]   0.15865358  0.78370189  0.02699149  0.00918335 363.32366541
[11,] 170.81748744  0.88369726 41.94086032 11.86782506   4.75586928
[12,]   3.88356742 16.23972661 10.50823459  2.27537972  10.65434611
[13,]   2.91243876  7.56993591  4.90401343  0.34693576  10.45467244
             [,11]       [,12]       [,13]
 [1,] 0.000000e+00    0.000000   0.0000000
 [2,] 0.000000e+00    0.000000   0.0000000
 [3,] 0.000000e+00    0.000000   0.0000000
 [4,] 0.000000e+00    0.000000   0.0000000
 [5,] 0.000000e+00    0.000000   0.0000000
 [6,] 1.708175e+02    3.883567   2.9124388
 [7,] 8.836973e-01   16.239727   7.5699359
 [8,] 4.194086e+01   10.508235   4.9040134
 [9,] 1.186783e+01    2.275380   0.3469358
[10,] 4.755869e+00   10.654346  10.4546724
[11,] 4.466207e+05 1356.317541  67.4930658
[12,] 1.356318e+03 2513.618033 116.5709449
[13,] 6.749307e+01  116.570945 353.1236854

 

 
************************** EXPECTED STANDARD ERRORS ************************
 
------------------------ Fixed Effects Parameters -------------------------
 
       Beta    StdError      RSE  
Kd     0.02  0.01298812 64.94061 %
CLD    0.20  0.08385527 41.92764 %
CLT    1.00  0.27308500 27.30850 %
RateT 40.00 12.95629788 32.39074 %
V      6.00  0.69627199 11.60453 %

 

 
------------------------- Variance of Inter-Subject Random Effects ----------------------
 
      Omega   StdError        RSE  
Kd     0.05 0.85621180 1712.42361 %
CLD    0.20 0.43672740  218.36370 %
CLT    0.25 0.16764869   67.05948 %
RateT  0.50 0.28318085   56.63617 %
V      0.10 0.05249098   52.49098 %

 
------------------------ Standard deviation of residual error ------------------------
 
           Sigma    StdError      RSE  
sig.interA  0.01 0.001531419 15.31419 %
sig.slopeA  0.20 0.020284041 10.14202 %
sig.slopeB  0.10 0.054737311 54.73731 %

 
******************************* DETERMINANT ********************************
 
1.894088e+23
 
******************************** CRITERION *********************************
 
61.74037
 

 

 
Time difference of 51.172 secs
sys.self 
    0.13 

 
