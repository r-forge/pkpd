PFIM 3.2  
 
Option: 1 

 
Project: Example Escalating Dose (Within ODE) Nonlinear Monkey PK - Dosing every 2nd day - only PK samples
 
Date: Thu Jan 13 11:31:20 2011
 

 
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
                                                                 Protocol
1 c=(0.0416666666666667, 1, 2, 2.04166666666667, 3, 4, 5, 6, 7, 8, 9, 10)
  subjects            condinit
1        1 c(0, RateT/CLT * V)

 
Sample times for response: B 
  Protocol subjects            condinit
1    c=(2)        1 c(0, RateT/CLT * V)

 
Total number of samples: 78
 
Associated criterion value: 33.3492
 
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

 

Now evaluating the Fisher Information Matrix for the 1 protocols generated 

 
**************************** OPTIMISED DESIGN *****************************
 

 
Optimised population design: 
Sample times for response: A 
                                                                   times freq
1 c(0.0416666666666667, 1, 2, 2.04166666666667, 3, 4, 5, 6, 7, 8, 9, 10)    1
  Subjects            condinit
1        6 c(0, RateT/CLT * V)

 
Sample times for response: B 
  times freq Subjects            condinit
1     2    1        6 c(0, RateT/CLT * V)

 

 
Associated optimised criterion: 33.3492
 

 
******************* POPULATION FISHER INFORMATION MATRIX ******************
 
             [,1]        [,2]         [,3]         [,4]         [,5]
 [1,] 272.4240034 27.00189382 -38.79976552 -0.544851855  1.128657618
 [2,]  27.0018938 65.28542134   5.89750992  0.021979188  1.596692914
 [3,] -38.7997655  5.89750992  19.88408879 -0.044115847 -0.035184656
 [4,]  -0.5448519  0.02197919  -0.04411585  0.006961562  0.001641010
 [5,]   1.1286576  1.59669291  -0.03518466  0.001641010  1.592583856
 [6,]   0.0000000  0.00000000   0.00000000  0.000000000  0.000000000
 [7,]   0.0000000  0.00000000   0.00000000  0.000000000  0.000000000
 [8,]   0.0000000  0.00000000   0.00000000  0.000000000  0.000000000
 [9,]   0.0000000  0.00000000   0.00000000  0.000000000  0.000000000
[10,]   0.0000000  0.00000000   0.00000000  0.000000000  0.000000000
[11,]   0.0000000  0.00000000   0.00000000  0.000000000  0.000000000
[12,]   0.0000000  0.00000000   0.00000000  0.000000000  0.000000000
[13,]   0.0000000  0.00000000   0.00000000  0.000000000  0.000000000
              [,6]        [,7]        [,8]         [,9]        [,10]
 [1,]   0.00000000 0.000000000  0.00000000  0.000000000   0.00000000
 [2,]   0.00000000 0.000000000  0.00000000  0.000000000   0.00000000
 [3,]   0.00000000 0.000000000  0.00000000  0.000000000   0.00000000
 [4,]   0.00000000 0.000000000  0.00000000  0.000000000   0.00000000
 [5,]   0.00000000 0.000000000  0.00000000  0.000000000   0.00000000
 [6,]   9.89531168 0.097213636  5.01807268  1.583272234   0.15286416
 [7,]   0.09721364 0.568291499  0.11593541  0.002576452   0.30593139
 [8,]   5.01807268 0.115935411 32.94808225  0.259494392   0.00371388
 [9,]   1.58327223 0.002576452  0.25949439 10.338848221   0.01292599
[10,]   0.15286416 0.305931391  0.00371388  0.012925993 273.92292056
[11,] 187.47532741 0.831367250 32.48778538  8.459648564   3.95694436
[12,]  14.57775028 7.652347064  9.59130523  2.107796790   6.58634954
[13,]   9.55358773 7.138653698  3.58490165  0.018055753  10.39153448
             [,11]       [,12]        [,13]
 [1,] 0.000000e+00    0.000000   0.00000000
 [2,] 0.000000e+00    0.000000   0.00000000
 [3,] 0.000000e+00    0.000000   0.00000000
 [4,] 0.000000e+00    0.000000   0.00000000
 [5,] 0.000000e+00    0.000000   0.00000000
 [6,] 1.874753e+02   14.577750   9.55358773
 [7,] 8.313673e-01    7.652347   7.13865370
 [8,] 3.248779e+01    9.591305   3.58490165
 [9,] 8.459649e+00    2.107797   0.01805575
[10,] 3.956944e+00    6.586350  10.39153448
[11,] 3.183335e+05 1947.346325  33.65556565
[12,] 1.947346e+03 1769.287397 109.46328827
[13,] 3.365557e+01  109.463288 271.49106519

 

 
************************** EXPECTED STANDARD ERRORS ************************
 
------------------------ Fixed Effects Parameters -------------------------
 
      Beta    StdError      RSE  
Kd     0.2  0.09482504 47.41252 %
CLD    0.2  0.14165623 70.82811 %
CLT    1.0  0.32018695 32.01870 %
RateT 40.0 15.09217327 37.73043 %
V      6.0  0.80280512 13.38009 %

 

 
------------------------- Variance of Inter-Subject Random Effects ----------------------
 
      Omega   StdError       RSE  
Kd     0.05 0.34560942 691.21883 %
CLD    0.20 1.65874754 829.37377 %
CLT    0.25 0.18159510  72.63804 %
RateT  0.50 0.31524038  63.04808 %
V      0.10 0.06046577  60.46577 %

 
------------------------ Standard deviation of residual error ------------------------
 
           Sigma    StdError      RSE  
sig.interA  0.01 0.001788128 17.88128 %
sig.slopeA  0.20 0.024698906 12.34945 %
sig.slopeB  0.10 0.075758893 75.75889 %

 
******************************* DETERMINANT ********************************
 
6.311212e+19
 
******************************** CRITERION *********************************
 
33.34921
 

 

 
Time difference of 8.641 secs
sys.self 
    0.08 

 
