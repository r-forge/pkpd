PFIM 3.2  
 
Option: 1 

 
Project: Example 3
 
Date: Mon Jan 10 21:39:23 2011
 

 
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

 
Initial Conditions at time 0 : 
 
0 0 

Error tolerance for solving differential equations system: RtolEQ = 1e-08 , AtolEQ = 1e-08 , Hmax =  0.5
 
Initial population design: 

 
Sample times for response: A 
                 Protocol subjects condinit
1 c=(0.5, 2, 30, 49, 180)      100  c(0, 0)

 
Sample times for response: B 
                  Protocol subjects condinit
1 c=(0.5, 2, 14, 110, 150)      100  c(0, 0)

 
Total number of samples: 1000
 
Associated criterion value: 1880.238
 
Identical sampling times for each response: FALSE
 
Random effect model: Trand =  2
 
Variance error model response A : ( 0 + 0.2 *f)^2
Variance error model response B : ( 0.1 + 0 *f)^2

 

Optimization step:  
 
Sampling windows for the response: A 
Window 1 : t= 0.5 2 10 30 49 158 180 
    Nb of sampling points to be taken in this window, n[ 1 ]= 5 
Maximum total number of points in one elementary protocol : 5 
Minimum total number of points in one elementary protocol : 5 

 
Sampling windows for the response: B 
Window 1 : t= 0.5 2 14 50 80 110 150 
    Nb of sampling points to be taken in this window, n[ 1 ]= 5 
Maximum total number of points in one elementary protocol : 5 
Minimum total number of points in one elementary protocol : 5 

 

Now evaluating the Fisher Information Matrix for the 441 protocols generated 

 
**************************** OPTIMISED DESIGN *****************************
 

 
Optimised population design: 
Sample times for response: A 
                   times freq Subjects condinit
1 c(0.5, 2, 30, 49, 180)    1      100  c(0, 0)

 
Sample times for response: B 
                   times freq Subjects condinit
1 c(0.5, 2, 50, 80, 110)    1      100  c(0, 0)

 

 
Associated optimised criterion: 1976.253
 

 
******************* POPULATION FISHER INFORMATION MATRIX ******************
 
             [,1]         [,2]          [,3]        [,4]       [,5]
 [1,]   2.5847770     5.813567     0.7773229   -11.39828   0.000000
 [2,]   5.8135667 57857.573248 -6821.4734747   355.31353   0.000000
 [3,]   0.7773229 -6821.473475   986.4310113  -152.86906   0.000000
 [4,] -11.3982828   355.313528  -152.8690604 37602.92368   0.000000
 [5,]   0.0000000     0.000000     0.0000000     0.00000 740.040476
 [6,]   0.0000000     0.000000     0.0000000     0.00000   0.169123
 [7,]   0.0000000     0.000000     0.0000000     0.00000   0.966871
 [8,]   0.0000000     0.000000     0.0000000     0.00000  69.564962
 [9,]   0.0000000     0.000000     0.0000000     0.00000   2.169749
              [,6]         [,7]        [,8]         [,9]
 [1,]   0.00000000   0.00000000     0.00000     0.000000
 [2,]   0.00000000   0.00000000     0.00000     0.000000
 [3,]   0.00000000   0.00000000     0.00000     0.000000
 [4,]   0.00000000   0.00000000     0.00000     0.000000
 [5,]   0.16912300   0.96687097    69.56496     2.169749
 [6,] 756.73852049   0.04244448    46.55538    12.482021
 [7,]   0.04244448 706.98993446    61.56690    97.162326
 [8,]  46.55538194  61.56689612 15135.76353   998.688975
 [9,]  12.48202146  97.16232571   998.68898 75637.159021

 

 
************************** EXPECTED STANDARD ERRORS ************************
 
------------------------ Fixed Effects Parameters -------------------------
 
       Beta    StdError       RSE  
V    12.200 0.623845226  5.113485 %
Vm    0.082 0.009701343 11.830906 %
km    0.370 0.074318588 20.086105 %
Alin  0.100 0.005164672  5.164672 %

 

 
------------------------- Variance of Inter-Subject Random Effects ----------------------
 
     Omega   StdError      RSE  
V     0.25 0.03676769 14.70708 %
Vm    0.25 0.03635537 14.54215 %
Alin  0.25 0.03761887 15.04755 %

 
------------------------ Standard deviation of residual error ------------------------
 
           Sigma    StdError      RSE  
sig.slopeA   0.2 0.008135709 4.067855 %
sig.interB   0.1 0.003637955 3.637955 %

 
******************************* DETERMINANT ********************************
 
4.598139e+29
 
******************************** CRITERION *********************************
 
1976.253
 

 

 
Time difference of 38.453 secs
sys.self 
    0.42 

 
