PFIM 3.2  
 
Option: 1 

 
Project: Example Escalating Dose (Within ODE) Nonlinear Monkey PK
 
Date: Wed Jan 12 14:08:37 2011
 

 
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
ndoses = 7
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
 
Associated criterion value: 54.8618
 
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

 

 
Associated optimised criterion: 54.8618
 

 
******************* POPULATION FISHER INFORMATION MATRIX ******************
 
              [,1]        [,2]        [,3]         [,4]        [,5]
 [1,] 2803.5256919  85.0449709 -18.7730636 -0.227960543 0.804984609
 [2,]   85.0444865 193.0161380 -16.2661317 -0.210658502 0.814262756
 [3,]  -18.7732564 -16.2660462  14.0787481 -0.121974140 0.314445566
 [4,]   -0.2279629  -0.2106575  -0.1219741  0.005988572 0.004029385
 [5,]    0.8049901   0.8142592   0.3144451  0.004029379 0.541878928
 [6,]    0.0000000   0.0000000   0.0000000  0.000000000 0.000000000
 [7,]    0.0000000   0.0000000   0.0000000  0.000000000 0.000000000
 [8,]    0.0000000   0.0000000   0.0000000  0.000000000 0.000000000
 [9,]    0.0000000   0.0000000   0.0000000  0.000000000 0.000000000
[10,]    0.0000000   0.0000000   0.0000000  0.000000000 0.000000000
[11,]    0.0000000   0.0000000   0.0000000  0.000000000 0.000000000
[12,]    0.0000000   0.0000000   0.0000000  0.000000000 0.000000000
[13,]    0.0000000   0.0000000   0.0000000  0.000000000 0.000000000
              [,6]        [,7]       [,8]       [,9]       [,10]        [,11]
 [1,] 0.000000e+00  0.00000000  0.0000000 0.00000000  0.00000000    0.0000000
 [2,] 0.000000e+00  0.00000000  0.0000000 0.00000000  0.00000000    0.0000000
 [3,] 0.000000e+00  0.00000000  0.0000000 0.00000000  0.00000000    0.0000000
 [4,] 0.000000e+00  0.00000000  0.0000000 0.00000000  0.00000000    0.0000000
 [5,] 0.000000e+00  0.00000000  0.0000000 0.00000000  0.00000000    0.0000000
 [6,] 1.047968e+03  0.96434745  1.1747718 0.27715496  0.07776056   19.2273200
 [7,] 9.643475e-01  4.96736394  0.8819522 0.23667618  0.07956252    7.8474815
 [8,] 1.174772e+00  0.88195217 16.5175956 1.98369217  0.29662764    6.4970387
 [9,] 2.771550e-01  0.23667618  1.9836922 7.65077181  0.07793242    1.7725722
[10,] 7.776056e-02  0.07956252  0.2966276 0.07793242 31.71233939    0.8200141
[11,] 1.922732e+01  7.84748154  6.4970387 1.77257217  0.82001411 2270.3981207
[12,] 5.990264e-01  0.71674966  1.3892505 0.09361394  0.11717336   18.5609103
[13,] 2.141977e+01 32.70704274 17.4969121 4.06075953  0.38670356  120.9466957
             [,12]        [,13]
 [1,] 0.000000e+00 0.000000e+00
 [2,] 0.000000e+00 0.000000e+00
 [3,] 0.000000e+00 0.000000e+00
 [4,] 0.000000e+00 0.000000e+00
 [5,] 0.000000e+00 0.000000e+00
 [6,] 5.990264e-01 2.141977e+01
 [7,] 7.167497e-01 3.270704e+01
 [8,] 1.389250e+00 1.749691e+01
 [9,] 9.361394e-02 4.060760e+00
[10,] 1.171734e-01 3.867036e-01
[11,] 1.856091e+01 1.209467e+02
[12,] 5.219127e+03 2.408911e+02
[13,] 2.408911e+02 1.229046e+04

 

 
************************** EXPECTED STANDARD ERRORS ************************
 
------------------------ Fixed Effects Parameters -------------------------
 
      Beta    StdError      RSE  
Kd     0.2  0.01911996  9.55998 %
CLD    0.2  0.08345070 41.72535 %
CLT    1.0  0.33788429 33.78843 %
RateT 40.0 15.74001577 39.35004 %
V      6.0  1.40514372 23.41906 %

 

 
------------------------- Variance of Inter-Subject Random Effects ----------------------
 
      Omega   StdError       RSE  
Kd     0.05 0.03089648  61.79296 %
CLD    0.20 0.45584840 227.92420 %
CLT    0.25 0.25122347 100.48939 %
RateT  0.50 0.36739335  73.47867 %
V      0.30 0.17759499  59.19833 %

 
------------------------ Standard deviation of residual error ------------------------
 
           Sigma    StdError       RSE  
sig.slopeA   0.2 0.021055557 10.527778 %
sig.slopeB   0.1 0.013848573 13.848573 %
sig.slopeC   0.1 0.009108264  9.108264 %

 
******************************* DETERMINANT ********************************
 
4.078638e+22
 
******************************** CRITERION *********************************
 
54.86184
 

 

 
Time difference of 21.015 secs
sys.self 
    0.06 

 
