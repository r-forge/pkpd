PFIM 3.2  
 
Option: 1 

 
Project: Example Escalating Dose (Within ODE) Nonlinear Monkey PK - Dosing every 2nd day - only PK samples
 
Date: Thu Jan 13 11:22:36 2011
 

 
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
1 c=(0.00694444444444444, 0.0416666666666667, 0.25, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  subjects            condinit
1        1 c(0, RateT/CLT * V)

 
Sample times for response: B 
  Protocol subjects            condinit
1    c=(2)        1 c(0, RateT/CLT * V)

 
Total number of samples: 90
 
Associated criterion value: 31.8147
 
Identical sampling times for each response: FALSE
 
Random effect model: Trand =  2
 
Variance error model response A : ( 0.01 + 0.2 *f)^2
Variance error model response B : ( 0 + 0.1 *f)^2

 

Optimization step:  
 
Sampling windows for the response: A 
Window 1 : t= 0.006944444 0.04166667 0.25 0.5 1 2 3 4 5 6 7 8 9 10 
    Nb of sampling points to be taken in this window, n[ 1 ]= 14 
Maximum total number of points in one elementary protocol : 14 
Minimum total number of points in one elementary protocol : 14 

 
Sampling windows for the response: B 
Window 1 : t= 2 
    Nb of sampling points to be taken in this window, n[ 1 ]= 1 
Maximum total number of points in one elementary protocol : 1 
Minimum total number of points in one elementary protocol : 1 

 

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
1     2    1        6 c(0, RateT/CLT * V)

 

 
Associated optimised criterion: 31.8147
 

 
******************* POPULATION FISHER INFORMATION MATRIX ******************
 
             [,1]        [,2]         [,3]         [,4]        [,5]
 [1,] 305.2419249 25.25657658 -38.98795580 -0.546195423 1.186154810
 [2,]  25.2565766 55.45717689   0.17061068 -0.044809464 2.132664301
 [3,] -38.9879558  0.17061068  16.56783585 -0.082769591 0.276178679
 [4,]  -0.5461954 -0.04480946  -0.08276959  0.006511042 0.005271246
 [5,]   1.1861548  2.13266430   0.27617868  0.005271246 1.563398350
 [6,]   0.0000000  0.00000000   0.00000000  0.000000000 0.000000000
 [7,]   0.0000000  0.00000000   0.00000000  0.000000000 0.000000000
 [8,]   0.0000000  0.00000000   0.00000000  0.000000000 0.000000000
 [9,]   0.0000000  0.00000000   0.00000000  0.000000000 0.000000000
[10,]   0.0000000  0.00000000   0.00000000  0.000000000 0.000000000
[11,]   0.0000000  0.00000000   0.00000000  0.000000000 0.000000000
[12,]   0.0000000  0.00000000   0.00000000  0.000000000 0.000000000
[13,]   0.0000000  0.00000000   0.00000000  0.000000000 0.000000000
              [,6]         [,7]         [,8]        [,9]       [,10]
 [1,]   0.00000000 0.000000e+00 0.000000e+00  0.00000000   0.0000000
 [2,]   0.00000000 0.000000e+00 0.000000e+00  0.00000000   0.0000000
 [3,]   0.00000000 0.000000e+00 0.000000e+00  0.00000000   0.0000000
 [4,]   0.00000000 0.000000e+00 0.000000e+00  0.00000000   0.0000000
 [5,]   0.00000000 0.000000e+00 0.000000e+00  0.00000000   0.0000000
 [6,]  12.42301770 8.505262e-02 5.066869e+00  1.59109035   0.1688356
 [7,]   0.08505262 4.100665e-01 9.702668e-05  0.01070874   0.5457908
 [8,]   5.06686899 9.702668e-05 2.287443e+01  0.91344068   0.2288240
 [9,]   1.59109035 1.070874e-02 9.134407e-01  9.04398305   0.1333730
[10,]   0.16883559 5.457908e-01 2.288240e-01  0.13337295 263.9751553
[11,] 303.00412168 3.288031e+00 7.744297e+01 18.08923072  14.5410686
[12,]  14.10797329 6.526717e+00 1.102907e+01  2.61882214  10.3240034
[13,]   9.23426346 5.876161e+00 6.844452e+00  0.22240276   7.9411332
             [,11]       [,12]       [,13]
 [1,] 0.000000e+00    0.000000   0.0000000
 [2,] 0.000000e+00    0.000000   0.0000000
 [3,] 0.000000e+00    0.000000   0.0000000
 [4,] 0.000000e+00    0.000000   0.0000000
 [5,] 0.000000e+00    0.000000   0.0000000
 [6,] 3.030041e+02   14.107973   9.2342635
 [7,] 3.288031e+00    6.526717   5.8761614
 [8,] 7.744297e+01   11.029069   6.8444516
 [9,] 1.808923e+01    2.618822   0.2224028
[10,] 1.454107e+01   10.324003   7.9411332
[11,] 6.297845e+05 2436.163213  66.2265234
[12,] 2.436163e+03 1595.065371 103.1969636
[13,] 6.622652e+01  103.196964 253.9664920

 

 
************************** EXPECTED STANDARD ERRORS ************************
 
------------------------ Fixed Effects Parameters -------------------------
 
      Beta   StdError      RSE  
Kd     0.2  0.0938023 46.90115 %
CLD    0.2  0.1417567 70.87834 %
CLT    1.0  0.3781025 37.81025 %
RateT 40.0 17.1520080 42.88002 %
V      6.0  0.8344141 13.90690 %

 

 
------------------------- Variance of Inter-Subject Random Effects ----------------------
 
      Omega   StdError       RSE  
Kd     0.05 0.30661649 613.23299 %
CLD    0.20 1.96161490 980.80745 %
CLT    0.25 0.21983733  87.93493 %
RateT  0.50 0.33653655  67.30731 %
V      0.10 0.06163413  61.63413 %

 
------------------------ Standard deviation of residual error ------------------------
 
           Sigma    StdError      RSE  
sig.interA  0.01 0.001270912 12.70912 %
sig.slopeA  0.20 0.026087677 13.04384 %
sig.slopeB  0.10 0.078132251 78.13225 %

 
******************************* DETERMINANT ********************************
 
3.421058e+19
 
******************************** CRITERION *********************************
 
31.81469
 

 

 
Time difference of 8.515 secs
sys.self 
    0.16 

 
