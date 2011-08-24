PFIM 3.2  
 
Option: 1 

 
Project: Example 5 
 
Date: Sun Jan 02 23:15:53 2011
 

 
**************************** INPUT SUMMARY ********************************
 
Analytical function model:  
 
dose/V * ka/(ka - (Cl/V)) * (exp(-(Cl/V) * t) - exp(-ka * t)) 
 

 
Initial population design: 

 
Sample times for response: A 
          Protocol subjects doses
1 c=(0.5, 2, 4, 8)       40    30

 
Total number of samples: 160
 
Associated criterion value: 1977.381
 
Identical sampling times for each response: TRUE
 
Number of occasions: 2
 
Random effect model: Trand =  2
 
Variance error model response A : ( 0.1 + 0 *f)^2

 

Optimization step:  
 
Sampling windows for the response: A 
Window 1 : t= 0.5 1 1.5 2 4 6 8 
    Nb of sampling points to be taken in this window, n[ 1 ]= 4 
Maximum total number of points in one elementary protocol : 4 
Minimum total number of points in one elementary protocol : 4 

 

Now evaluating the Fisher Information Matrix for the 35 protocols generated 

 
**************************** OPTIMISED DESIGN *****************************
 

 
Optimised population design: 
Sample times for response: A 
            times freq Subjects doses
1 c(0.5, 2, 6, 8)    1       40    30

 

 
Associated optimised criterion: 2068.720
 

 
******************* POPULATION FISHER INFORMATION MATRIX ******************
 
            [,1]       [,2]      [,3]         [,4]         [,5]         [,6]
 [1,] 342.070959 -11.909950  2.294299    0.0000000    0.0000000    0.0000000
 [2,] -11.909950  29.371134  0.624963    0.0000000    0.0000000    0.0000000
 [3,]   2.294299   0.624963 98.030214    0.0000000    0.0000000    0.0000000
 [4,]   0.000000   0.000000  0.000000 1462.6567613   21.7203074    0.2631903
 [5,]   0.000000   0.000000  0.000000   21.7203074 1618.1679873    0.2392295
 [6,]   0.000000   0.000000  0.000000    0.2631903    0.2392295 1921.9845543
 [7,]   0.000000   0.000000  0.000000  731.3283806   10.8601537    0.1315952
 [8,]   0.000000   0.000000  0.000000   10.8601537  809.0839937    0.1196147
 [9,]   0.000000   0.000000  0.000000    0.1315952    0.1196147  960.9922772
[10,]   0.000000   0.000000  0.000000  828.6261286  553.4112361   56.3448400
              [,7]         [,8]         [,9]       [,10]
 [1,] 0.000000e+00 0.000000e+00 0.000000e+00     0.00000
 [2,] 0.000000e+00 0.000000e+00 0.000000e+00     0.00000
 [3,] 0.000000e+00 0.000000e+00 0.000000e+00     0.00000
 [4,] 7.313284e+02 1.086015e+01 1.315952e-01   828.62613
 [5,] 1.086015e+01 8.090840e+02 1.196147e-01   553.41124
 [6,] 1.315952e-01 1.196147e-01 9.609923e+02    56.34484
 [7,] 1.673941e+04 4.697134e+03 4.315501e+01  2608.00290
 [8,] 4.697134e+03 2.519417e+04 6.569500e+01  1889.93103
 [9,] 4.315501e+01 6.569500e+01 4.381684e+04   926.11512
[10,] 2.608003e+03 1.889931e+03 9.261151e+02 20551.58326

 

 
************************** EXPECTED STANDARD ERRORS ************************
 
------------------------ Fixed Effects Parameters -------------------------
 
   Beta   StdError      RSE  
ka  1.0 0.05445931 5.445931 %
V   3.5 0.18585115 5.310033 %
Cl  2.0 0.10101646 5.050823 %

 

 
------------------------- Variance of Inter-Subject Random Effects ----------------------
 
   Omega   StdError      RSE  
ka  0.09 0.02669810 29.66456 %
V   0.09 0.02517532 27.97258 %
Cl  0.09 0.02293649 25.48499 %

 
------------------------- Variance of Inter-Occasion Random Effects ----------------------
 
    Gamma    StdError      RSE  
ka 0.0225 0.008083050 35.92466 %
V  0.0225 0.006537789 29.05684 %
Cl 0.0225 0.004805936 21.35972 %

 
------------------------ Standard deviation of residual error ------------------------
 
           Sigma    StdError      RSE  
sig.interA   0.1 0.007152717 7.152717 %

 
******************************* DETERMINANT ********************************
 
1.435548e+33
 
******************************** CRITERION *********************************
 
2068.720
 

 

 
Time difference of 0.937 secs
sys.self 
    0.01 

 
