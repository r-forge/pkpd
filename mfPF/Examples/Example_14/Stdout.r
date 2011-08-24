PFIM 3.2  
 
Option: 1 

 
Project: Example Escalating Dose (Within ODE) Nonlinear Monkey PK - Dosing every 2nd day
 
Date: Sun Jan 23 23:56:38 2011
 

 
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

Error tolerance for solving differential equations system: RtolEQ = 1e-08 , AtolEQ = 1e-08 , Hmax =  0.01
 
Initial population design: 

 
Sample times for response: A 
                                                                               Protocol
1 c=(0.00694444444444444, 0.0416666666666667, 0.25, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  subjects            condinit
1        6 c(0, RateT/CLT * V)

 
Sample times for response: B 
            Protocol subjects            condinit
1 c=(1, 3, 5, 7, 10)        6 c(0, RateT/CLT * V)

 
Sample times for response: C 
                                                          Protocol subjects
1 c=(0.0416666666666667, 0.25, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)        6
             condinit
1 c(0, RateT/CLT * V)

 
Total number of samples: 192
 
Associated criterion value: 54.5027
 
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

 

 
Associated optimised criterion: 54.5027
 

 
******************* POPULATION FISHER INFORMATION MATRIX ******************
 
              [,1]         [,2]         [,3]         [,4]        [,5]
 [1,] 2781.5880909  86.25705201 -11.17965997 -0.134020856 0.739413765
 [2,]   86.2570519 105.83268269  -0.49415579 -0.013176094 0.425567408
 [3,]  -11.1796740  -0.49415510  19.61118670 -0.053382979 0.164174659
 [4,]   -0.1340210  -0.01317608  -0.05338298  0.006839744 0.002151662
 [5,]    0.7394142   0.42556739   0.16417462  0.002151662 0.545332125
 [6,]    0.0000000   0.00000000   0.00000000  0.000000000 0.000000000
 [7,]    0.0000000   0.00000000   0.00000000  0.000000000 0.000000000
 [8,]    0.0000000   0.00000000   0.00000000  0.000000000 0.000000000
 [9,]    0.0000000   0.00000000   0.00000000  0.000000000 0.000000000
[10,]    0.0000000   0.00000000   0.00000000  0.000000000 0.000000000
[11,]    0.0000000   0.00000000   0.00000000  0.000000000 0.000000000
[12,]    0.0000000   0.00000000   0.00000000  0.000000000 0.000000000
[13,]    0.0000000   0.00000000   0.00000000  0.000000000 0.000000000
              [,6]         [,7]         [,8]         [,9]       [,10]
 [1,] 0.000000e+00 0.000000e+00 0.000000e+00 0.0000000000  0.00000000
 [2,] 0.000000e+00 0.000000e+00 0.000000e+00 0.0000000000  0.00000000
 [3,] 0.000000e+00 0.000000e+00 0.000000e+00 0.0000000000  0.00000000
 [4,] 0.000000e+00 0.000000e+00 0.000000e+00 0.0000000000  0.00000000
 [5,] 0.000000e+00 0.000000e+00 0.000000e+00 0.0000000000  0.00000000
 [6,] 1.031631e+03 9.920372e-01 4.166165e-01 0.0957952701  0.06560797
 [7,] 9.920372e-01 1.493408e+00 8.139653e-04 0.0009259164  0.02173291
 [8,] 4.166165e-01 8.139653e-04 3.204989e+01 0.3799656666  0.08085994
 [9,] 9.579527e-02 9.259164e-04 3.799657e-01 9.9801822989  0.02222231
[10,] 6.560797e-02 2.173291e-02 8.085994e-02 0.0222223125 32.11780970
[11,] 2.337026e+01 5.772247e+00 4.818989e+00 1.2659603069  0.80625481
[12,] 4.369922e-01 6.531740e-01 7.273602e-01 0.0263707258  0.11875234
[13,] 2.642391e+01 2.302601e+01 2.079241e+01 4.9431363738  1.10238247
             [,11]        [,12]        [,13]
 [1,]    0.0000000 0.000000e+00     0.000000
 [2,]    0.0000000 0.000000e+00     0.000000
 [3,]    0.0000000 0.000000e+00     0.000000
 [4,]    0.0000000 0.000000e+00     0.000000
 [5,]    0.0000000 0.000000e+00     0.000000
 [6,]   23.3702581 4.369922e-01    26.423912
 [7,]    5.7722467 6.531740e-01    23.026012
 [8,]    4.8189892 7.273602e-01    20.792406
 [9,]    1.2659603 2.637073e-02     4.943136
[10,]    0.8062548 1.187523e-01     1.102382
[11,] 1931.7550601 1.452654e+01   164.099417
[12,]   14.5265386 5.259639e+03   246.320633
[13,]  164.0994169 2.463206e+02 12009.091217

 

 
************************** EXPECTED STANDARD ERRORS ************************
 
------------------------ Fixed Effects Parameters -------------------------
 
      Beta    StdError      RSE  
Kd     0.2  0.01924026  9.62013 %
CLD    0.2  0.09861087 49.30543 %
CLT    1.0  0.22895189 22.89519 %
RateT 40.0 12.24390195 30.60975 %
V      6.0  1.35950927 22.65849 %

 

 
------------------------- Variance of Inter-Subject Random Effects ----------------------
 
      Omega   StdError       RSE  
Kd     0.05 0.03114738  62.29476 %
CLD    0.20 0.83537649 417.68824 %
CLT    0.25 0.17681137  70.72455 %
RateT  0.50 0.31665477  63.33095 %
V      0.30 0.17645458  58.81819 %

 
------------------------ Standard deviation of residual error ------------------------
 
           Sigma    StdError       RSE  
sig.slopeA   0.2 0.022895011 11.447506 %
sig.slopeB   0.1 0.013795396 13.795396 %
sig.slopeC   0.1 0.009274746  9.274746 %

 
******************************* DETERMINANT ********************************
 
3.744825e+22
 
******************************** CRITERION *********************************
 
54.50267
 

 

 
Time difference of 1.3052 mins
sys.self 
    0.11 

 
