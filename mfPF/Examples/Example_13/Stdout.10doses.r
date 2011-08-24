PFIM 3.2  
 
Option: 1 

 
Project: Example Escalating Dose (Within ODE) Nonlinear Monkey PK
 
Date: Wed Jan 12 14:06:19 2011
 

 
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
ndoses = 10
startt= 0                                        #days start time for dosing
tau   = 1                                        #days interdose interval of antibody

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
1 c=(0.00694444444444444, 0.0416666666666667, 0.25, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  subjects            condinit
1        1 c(0, RateT/CLT * V)

 
Sample times for response: B 
            Protocol subjects            condinit
1 c=(1, 3, 5, 7, 10)        1 c(0, RateT/CLT * V)

 
Sample times for response: C 
                                                          Protocol subjects
1 c=(0.0416666666666667, 0.25, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)        1
             condinit
1 c(0, RateT/CLT * V)

 
Total number of samples: 192
 
Associated criterion value: 36.9486
 
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
  freq Subjects            condinit
1    1        6 c(0, RateT/CLT * V)

 
Sample times for response: B 
              times freq Subjects            condinit
1 c(1, 3, 5, 7, 10)    1        6 c(0, RateT/CLT * V)

 
Sample times for response: C 
                                                            times freq Subjects
1 c(0.0416666666666667, 0.25, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)    1        6
             condinit
1 c(0, RateT/CLT * V)

 

 
Associated optimised criterion: 36.9486
 

 
******************* POPULATION FISHER INFORMATION MATRIX ******************
 
              [,1]        [,2]        [,3]         [,4]        [,5]
 [1,] 2800.9806092  66.9327129 -19.6668110 -0.239147885 0.839855528
 [2,]   66.9321802  68.9090336 -20.2099981 -0.259950957 0.976491717
 [3,]  -19.6669741 -20.2099297  15.1018011 -0.109127144 0.278885208
 [4,]   -0.2391499  -0.2599501  -0.1091271  0.006149898 0.003582800
 [5,]    0.8398600   0.9764886   0.2788848  0.003582795 0.543110410
 [6,]    0.0000000   0.0000000   0.0000000  0.000000000 0.000000000
 [7,]    0.0000000   0.0000000   0.0000000  0.000000000 0.000000000
 [8,]    0.0000000   0.0000000   0.0000000  0.000000000 0.000000000
 [9,]    0.0000000   0.0000000   0.0000000  0.000000000 0.000000000
[10,]    0.0000000   0.0000000   0.0000000  0.000000000 0.000000000
[11,]    0.0000000   0.0000000   0.0000000  0.000000000 0.000000000
[12,]    0.0000000   0.0000000   0.0000000  0.000000000 0.000000000
[13,]    0.0000000   0.0000000   0.0000000  0.000000000 0.000000000
              [,6]      [,7]       [,8]       [,9]       [,10]        [,11]
 [1,] 0.000000e+00 0.0000000  0.0000000 0.00000000  0.00000000    0.0000000
 [2,] 0.000000e+00 0.0000000  0.0000000 0.00000000  0.00000000    0.0000000
 [3,] 0.000000e+00 0.0000000  0.0000000 0.00000000  0.00000000    0.0000000
 [4,] 0.000000e+00 0.0000000  0.0000000 0.00000000  0.00000000    0.0000000
 [5,] 0.000000e+00 0.0000000  0.0000000 0.00000000  0.00000000    0.0000000
 [6,] 1.046066e+03 0.5973270  1.2892889 0.30502504  0.08464333   19.9233274
 [7,] 5.973270e-01 0.6331273  1.3614755 0.36039617  0.11442396    2.7919109
 [8,] 1.289289e+00 1.3614755 19.0053664 1.58783120  0.23333056    6.9468255
 [9,] 3.050250e-01 0.3603962  1.5878312 8.06853116  0.06161489    1.8877008
[10,] 8.464333e-02 0.1144240  0.2333306 0.06161489 31.85664308    0.8611326
[11,] 1.992333e+01 2.7919109  6.9468255 1.88770083  0.86113263 2282.0680295
[12,] 5.962592e-01 0.7156016  1.4704760 0.10433732  0.11229064   18.5417645
[13,] 2.148974e+01 7.0364700 16.6135507 3.77785147  0.22019069  120.0525340
             [,12]        [,13]
 [1,]    0.0000000 0.000000e+00
 [2,]    0.0000000 0.000000e+00
 [3,]    0.0000000 0.000000e+00
 [4,]    0.0000000 0.000000e+00
 [5,]    0.0000000 0.000000e+00
 [6,]    0.5962592 2.148974e+01
 [7,]    0.7156016 7.036470e+00
 [8,]    1.4704760 1.661355e+01
 [9,]    0.1043373 3.777851e+00
[10,]    0.1122906 2.201907e-01
[11,]   18.5417645 1.200525e+02
[12,] 5219.2609918 2.404320e+02
[13,]  240.4320335 1.251031e+04

 

 
************************** EXPECTED STANDARD ERRORS ************************
 
------------------------ Fixed Effects Parameters -------------------------
 
      Beta    StdError        RSE  
Kd     0.2  0.01911952   9.559758 %
CLD    0.2  0.43884246 219.421231 %
CLT    1.0  0.90785914  90.785914 %
RateT 40.0 37.32930141  93.323254 %
V      6.0  2.00328855  33.388143 %

 

 
------------------------- Variance of Inter-Subject Random Effects ----------------------
 
      Omega   StdError       RSE  
Kd     0.05 0.03092917  61.85835 %
CLD    0.20 1.38330584 691.65292 %
CLT    0.25 0.25005405 100.02162 %
RateT  0.50 0.35754820  71.50964 %
V      0.30 0.17723198  59.07733 %

 
------------------------ Standard deviation of residual error ------------------------
 
           Sigma    StdError       RSE  
sig.slopeA   0.2 0.020994630 10.497315 %
sig.slopeB   0.1 0.013848862 13.848862 %
sig.slopeC   0.1 0.008973776  8.973776 %

 
******************************* DETERMINANT ********************************
 
2.392039e+20
 
******************************** CRITERION *********************************
 
36.94856
 

 

 
Time difference of 24.094 secs
sys.self 
    0.14 

 
