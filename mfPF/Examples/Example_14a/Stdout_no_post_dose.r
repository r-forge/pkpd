PFIM 3.2  
 
Option: 1 

 
Project: Example Escalating Dose (Within ODE) Nonlinear Monkey PK - Dosing every 2nd day - only PK samples
 
Date: Thu Jan 13 11:27:11 2011
 

 
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
StartDose = 0.06 # For multi-dosing...
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
return(list(c(dTD,dTT),c(CFD,CTT)))
}

 
Initial Conditions at time -2 : 
 
0 RateT/CLT * V 

Error tolerance for solving differential equations system: RtolEQ = 1e-08 , AtolEQ = 1e-08 , Hmax =  0.01
 
Initial population design: 

 
Sample times for response: A 
                           Protocol subjects            condinit
1 c=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)        1 c(0, RateT/CLT * V)

 
Sample times for response: B 
  Protocol subjects            condinit
1    c=(2)        1 c(0, RateT/CLT * V)

 
Total number of samples: 66
 
Associated criterion value: 26.6184
 
Identical sampling times for each response: FALSE
 
Random effect model: Trand =  2
 
Variance error model response A : ( 0.01 + 0.2 *f)^2
Variance error model response B : ( 0 + 0.1 *f)^2

 

Optimization step:  
 
Sampling windows for the response: A 
Window 1 : t= 1 2 3 4 5 6 7 8 9 10 
    Nb of sampling points to be taken in this window, n[ 1 ]= 10 
Maximum total number of points in one elementary protocol : 10 
Minimum total number of points in one elementary protocol : 10 

 
Sampling windows for the response: B 
Window 1 : t= 2 
    Nb of sampling points to be taken in this window, n[ 1 ]= 1 
Maximum total number of points in one elementary protocol : 1 
Minimum total number of points in one elementary protocol : 1 

 

Now evaluating the Fisher Information Matrix for the 1 protocols generated 

 
**************************** OPTIMISED DESIGN *****************************
 

 
Optimised population design: 
Sample times for response: A 
                             times freq Subjects            condinit
1 c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)    1        6 c(0, RateT/CLT * V)

 
Sample times for response: B 
  times freq Subjects            condinit
1     2    1        6 c(0, RateT/CLT * V)

 

 
Associated optimised criterion: 26.6184
 

 
******************* POPULATION FISHER INFORMATION MATRIX ******************
 
             [,1]        [,2]         [,3]         [,4]        [,5]        [,6]
 [1,] 256.5296966 25.06427781 -40.34603017 -0.563282797 1.253331221   0.0000000
 [2,]  25.0642778 55.45629309   0.16506987 -0.044879012 2.132941612   0.0000000
 [3,] -40.3460302  0.16506987  16.52970753 -0.083249078 0.278069443   0.0000000
 [4,]  -0.5632828 -0.04487901  -0.08324908  0.006505012 0.005295019   0.0000000
 [5,]   1.2533312  2.13294161   0.27806944  0.005295019 1.563304503   0.0000000
 [6,]   0.0000000  0.00000000   0.00000000  0.000000000 0.000000000   8.7743314
 [7,]   0.0000000  0.00000000   0.00000000  0.000000000 0.000000000   0.0837624
 [8,]   0.0000000  0.00000000   0.00000000  0.000000000 0.000000000   5.4260072
 [9,]   0.0000000  0.00000000   0.00000000  0.000000000 0.000000000   1.6922000
[10,]   0.0000000  0.00000000   0.00000000  0.000000000 0.000000000   0.1885007
[11,]   0.0000000  0.00000000   0.00000000  0.000000000 0.000000000 127.0404370
[12,]   0.0000000  0.00000000   0.00000000  0.000000000 0.000000000  13.5500717
[13,]   0.0000000  0.00000000   0.00000000  0.000000000 0.000000000   9.2692729
              [,7]         [,8]        [,9]       [,10]        [,11]
 [1,] 0.000000e+00 0.000000e+00  0.00000000   0.0000000 0.000000e+00
 [2,] 0.000000e+00 0.000000e+00  0.00000000   0.0000000 0.000000e+00
 [3,] 0.000000e+00 0.000000e+00  0.00000000   0.0000000 0.000000e+00
 [4,] 0.000000e+00 0.000000e+00  0.00000000   0.0000000 0.000000e+00
 [5,] 0.000000e+00 0.000000e+00  0.00000000   0.0000000 0.000000e+00
 [6,] 8.376240e-02 5.426007e+00  1.69220005   0.1885007 1.270404e+02
 [7,] 4.100534e-01 9.082687e-05  0.01074200   0.5459328 3.285896e+00
 [8,] 9.082687e-05 2.276927e+01  0.92405454   0.2319678 7.397616e+01
 [9,] 1.074200e-02 9.240545e-01  9.02723896   0.1345787 1.720721e+01
[10,] 5.459328e-01 2.319678e-01  0.13457868 263.9434647 1.424104e+01
[11,] 3.285896e+00 7.397616e+01 17.20721362  14.2410399 2.139508e+05
[12,] 6.526989e+00 1.104085e+01  2.62021866  10.3291930 9.935293e+02
[13,] 5.876233e+00 6.840556e+00  0.22204745   7.9423437 6.618101e+01
            [,12]       [,13]
 [1,]    0.000000   0.0000000
 [2,]    0.000000   0.0000000
 [3,]    0.000000   0.0000000
 [4,]    0.000000   0.0000000
 [5,]    0.000000   0.0000000
 [6,]   13.550072   9.2692729
 [7,]    6.526989   5.8762335
 [8,]   11.040853   6.8405562
 [9,]    2.620219   0.2220474
[10,]   10.329193   7.9423437
[11,]  993.529327  66.1810100
[12,] 1589.828896 103.1917533
[13,]  103.191753 253.9632323

 

 
************************** EXPECTED STANDARD ERRORS ************************
 
------------------------ Fixed Effects Parameters -------------------------
 
      Beta   StdError      RSE  
Kd     0.2  0.1391814 69.59068 %
CLD    0.2  0.1450389 72.51943 %
CLT    1.0  0.5008903 50.08903 %
RateT 40.0 21.6861778 54.21544 %
V      6.0  0.8474570 14.12428 %

 

 
------------------------- Variance of Inter-Subject Random Effects ----------------------
 
      Omega   StdError       RSE  
Kd     0.05 0.38149055 762.98111 %
CLD    0.20 1.96431614 982.15807 %
CLT    0.25 0.22717517  90.87007 %
RateT  0.50 0.33930093  67.86019 %
V      0.10 0.06163792  61.63792 %

 
------------------------ Standard deviation of residual error ------------------------
 
           Sigma    StdError      RSE  
sig.interA  0.01 0.002173729 21.73729 %
sig.slopeA  0.20 0.026122191 13.06110 %
sig.slopeB  0.10 0.078604479 78.60448 %

 
******************************* DETERMINANT ********************************
 
3.367893e+18
 
******************************** CRITERION *********************************
 
26.61837
 

 

 
Time difference of 8.562 secs
sys.self 
    0.03 

 
