PFIM 3.2  
 
Option: 1 

 
Project: Example Single Dose Nonlinear Monkey PK
 
Date: Wed Jan 12 13:53:18 2011
 

 
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
input = 0

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

 
Initial Conditions at time 0 : 
 
5 * 70/0.15 RateT/CLT * V 

Error tolerance for solving differential equations system: RtolEQ = 1e-08 , AtolEQ = 1e-08 , Hmax =  0.5
 
Initial population design: 

 
Sample times for response: A 
                                                                                     Protocol
1 c=(0.00694444444444444, 0.0416666666666667, 0.25, 0.5, 2, 4, 7, 14, 21, 35, 49, 63, 77, 91)
  subjects                      condinit
1        1 c(5 * 70/0.15, RateT/CLT * V)

 
Sample times for response: B 
                                Protocol subjects                      condinit
1 c=(0.00694444444444444, 7, 35, 77, 91)        1 c(5 * 70/0.15, RateT/CLT * V)

 
Sample times for response: C 
                                      Protocol subjects
1 c=(0.5, 2, 4, 7, 14, 21, 35, 49, 63, 77, 91)        1
                       condinit
1 c(5 * 70/0.15, RateT/CLT * V)

 
Total number of samples: 180
 
Associated criterion value: 86.477
 
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
  freq Subjects                      condinit
1    1        6 c(5 * 70/0.15, RateT/CLT * V)

 
Sample times for response: B 
                                  times freq Subjects
1 c(0.00694444444444444, 7, 35, 77, 91)    1        6
                       condinit
1 c(5 * 70/0.15, RateT/CLT * V)

 
Sample times for response: C 
                                        times freq Subjects
1 c(0.5, 2, 4, 7, 14, 21, 35, 49, 63, 77, 91)    1        6
                       condinit
1 c(5 * 70/0.15, RateT/CLT * V)

 

 
Associated optimised criterion: 86.477
 

 
******************* POPULATION FISHER INFORMATION MATRIX ******************
 
               [,1]        [,2]        [,3]          [,4]          [,5]
 [1,] 2627.13276412  55.4126931 -4.00795388 -0.0164942068  0.9963613429
 [2,]   55.41269330 653.5176989  8.43770536  0.1052993052  0.1176559126
 [3,]   -4.00795390   8.4377054 22.99224187 -0.0110182978 -0.0360946983
 [4,]   -0.01649421   0.1052993 -0.01101830  0.0073654493 -0.0004338691
 [5,]    0.99636134   0.1176559 -0.03609470 -0.0004338691  0.5486917298
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
 [6,] 9.202435e+02 4.094089e-01  0.053545648 1.450981e-03 1.191283e-01
 [7,] 4.094089e-01 5.694472e+01  0.237316240 5.913570e-02 1.661150e-03
 [8,] 5.354565e-02 2.373162e-01 44.053598846 1.618705e-02 3.908482e-03
 [9,] 1.450981e-03 5.913570e-02  0.016187051 1.157330e+01 9.035634e-04
[10,] 1.191283e-01 1.661150e-03  0.003908482 9.035634e-04 3.251476e+01
[11,] 3.650744e+01 4.434568e+00  1.247285921 2.882486e-01 1.001322e+00
[12,] 1.505747e+01 3.945074e-04  0.830742331 5.386651e-03 3.597562e-03
[13,] 3.990257e+01 2.256136e+01  5.140978842 1.207955e+00 2.798562e-01
             [,11]        [,12]        [,13]
 [1,]    0.0000000 0.000000e+00    0.0000000
 [2,]    0.0000000 0.000000e+00    0.0000000
 [3,]    0.0000000 0.000000e+00    0.0000000
 [4,]    0.0000000 0.000000e+00    0.0000000
 [5,]    0.0000000 0.000000e+00    0.0000000
 [6,]   36.5074408 1.505747e+01   39.9025660
 [7,]    4.4345682 3.945074e-04   22.5613641
 [8,]    1.2472859 8.307423e-01    5.1409788
 [9,]    0.2882486 5.386651e-03    1.2079550
[10,]    1.0013220 3.597562e-03    0.2798562
[11,] 2358.6828930 3.513334e+00  224.1122486
[12,]    3.5133339 5.041008e+03  270.4592594
[13,]  224.1122486 2.704593e+02 9335.4555481

 

 
************************** EXPECTED STANDARD ERRORS ************************
 
------------------------ Fixed Effects Parameters -------------------------
 
      Beta    StdError      RSE  
Kd     0.2  0.01953806  9.76903 %
CLD    0.2  0.03929801 19.64901 %
CLT    1.0  0.20919153 20.91915 %
RateT 40.0 11.67130665 29.17827 %
V      6.0  1.35059410 22.50990 %

 

 
------------------------- Variance of Inter-Subject Random Effects ----------------------
 
      Omega   StdError      RSE  
Kd     0.05 0.03297806 65.95612 %
CLD    0.20 0.13259032 66.29516 %
CLT    0.25 0.15067138 60.26855 %
RateT  0.50 0.29395171 58.79034 %
V      0.30 0.17537297 58.45766 %

 
------------------------ Standard deviation of residual error ------------------------
 
           Sigma   StdError      RSE  
sig.slopeA   0.2 0.02062141 10.31071 %
sig.slopeB   0.1 0.01409578 14.09578 %
sig.slopeC   0.1 0.01037563 10.37563 %

 
******************************* DETERMINANT ********************************
 
1.512541e+25
 
******************************** CRITERION *********************************
 
86.47699
 

 

 
Time difference of 4.156 secs
sys.self 
    0.13 

 
