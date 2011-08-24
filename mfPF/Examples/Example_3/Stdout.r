PFIM 3.2 Option 1 
 
Project:  Example 3
 
Date:  Tue Jan 11 18:39:33 2011
 

 
**************************** INPUT SUMMARY ********************************
 
Differential Equations form of the model:  
 
function(t,y,p){
V<-p[1]
Vm<-p[2]
km<-p[3]
Alin<-p[4]
pk<-y[1:1]
pd<-y[2:2]
conc<-y[1]
if(t<=1){
dpk1<-(100/(1*V))+(-Vm)*pk[1]/(km*V+pk[1])}
else{
dpk1<-(-Vm)*pk[1]/(km*V+pk[1])}
dpd1<-0
pdIm<-Alin*conc
return(list(c(dpk1,dpd1),c(pk[1],pdIm)))
}

 
Population design:  
Sample times for response: A 
                   times subjects
1 c(0.5, 2, 30, 49, 180)      100

 
Sample times for response: B 
                    times subjects
1 c(0.5, 2, 14, 110, 150)      100

 
Variance error model response A : ( 0 + 0.2 *f)^2
Variance error model response B : ( 0.1 + 0 *f)^2

 
Initial Conditions at time 0 : 
 
0 0 

 
Random effect model: Trand =  2
 
Error tolerance for solving differential equations system: RtolEQ = 1e-08 , AtolEQ = 1e-08 , Hmax =  0.5
 

 
Computation of the Fisher information matrix: option =  1
 
******************* POPULATION FISHER INFORMATION MATRIX ******************
 
             [,1]         [,2]          [,3]        [,4]        [,5]         [,6]
 [1,]   2.5865260     5.374266     0.8128864   -11.44898   0.0000000   0.00000000
 [2,]   5.3742659 57945.602398 -6866.1966934   435.14502   0.0000000   0.00000000
 [3,]   0.8128864 -6866.196693   930.6251972   -26.96431   0.0000000   0.00000000
 [4,] -11.4489834   435.145019   -26.9643093 37662.89817   0.0000000   0.00000000
 [5,]   0.0000000     0.000000     0.0000000     0.00000 741.0423198   0.14452922
 [6,]   0.0000000     0.000000     0.0000000     0.00000   0.1445292 759.04299755
 [7,]   0.0000000     0.000000     0.0000000     0.00000   0.9754916   0.06365987
 [8,]   0.0000000     0.000000     0.0000000     0.00000  67.4492295  41.68128936
 [9,]   0.0000000     0.000000     0.0000000     0.00000   4.0750783  16.64379328
              [,7]        [,8]         [,9]
 [1,]   0.00000000     0.00000     0.000000
 [2,]   0.00000000     0.00000     0.000000
 [3,]   0.00000000     0.00000     0.000000
 [4,]   0.00000000     0.00000     0.000000
 [5,]   0.97549156    67.44923     4.075078
 [6,]   0.06365987    41.68129    16.643793
 [7,] 709.24694938    60.48191    93.895495
 [8,]  60.48190935 15215.85874  1305.080765
 [9,]  93.89549451  1305.08077 74085.359061

 
************************** EXPECTED STANDARD ERRORS ************************
 
------------------------ Fixed Effects Parameters -------------------------
 
       Beta    StdError       RSE  
V    12.200 0.624466721  5.118580 %
Vm    0.082 0.011757152 14.337990 %
km    0.370 0.092773757 25.073988 %
Alin  0.100 0.005157086  5.157086 %

 
------------------------- Variance of Inter-Subject Random Effects ----------------------
 
     Omega   StdError      RSE  
V     0.25 0.03674230 14.69692 %
Vm    0.25 0.03629946 14.51978 %
Alin  0.25 0.03755845 15.02338 %

 
------------------------ Standard deviation of residual error ------------------------ 
 
           Sigma    StdError      RSE  
sig.slopeA   0.2 0.008116514 4.058257 %
sig.interB   0.1 0.003677015 3.677015 %

 
******************************* DETERMINANT ********************************
 
2.937072e+29
 
******************************** CRITERION *********************************
 
1880.238
 

 
Time difference of 9.797 secs
