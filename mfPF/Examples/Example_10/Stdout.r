PFIM 3.2  
 
Option: 1 

 
Project: Example 10
 
Date: Sun Jan 02 23:58:41 2011
 

 
**************************** INPUT SUMMARY ********************************
 
Analytical function model:  
 
dose/V * ka/(ka - (Cl/V)) * (exp(-(Cl/V) * t) - exp(-ka * t)) 
 

 
Initial population design: 

 
Sample times for response: A 
             Protocol subjects doses
1 c=(0.5, 2, 4, 6, 8)       20    30
2   c=(1, 2, 3, 4, 5)       20    50

 
Total number of samples: 200
 
Associated criterion value: 555.091
 
Identical sampling times for each response: FALSE
 
Number of occasions: 2
 
Random effect model: Trand =  2
 
Variance error model response A : ( 0.5 + 0.1 *f)^2

 
Covariate model :  
 
	NB: Covariates are additive on log parameters
 
		Covariate 1 : Sex ( V ) 
    Categories References Proportions
(1)          M          *         0.5
(2)          F                    0.5


Optimization step:  
 
Sampling windows for the response: A 
Window 1 : t= 0.5 1 1.5 2 3 4 5 6 8 
    Nb of sampling points to be taken in this window, n[ 1 ]= 5 
Maximum total number of points in one elementary protocol : 5 
Minimum total number of points in one elementary protocol : 5 

 

Now evaluating the Fisher Information Matrix for the 252 protocols generated 

 
**************************** OPTIMISED DESIGN *****************************
 

 
Optimised population design: 
Sample times for response: A 
                 times        freq    Subjects doses
1 c(0.5, 1.5, 2, 6, 8) 0.492355997 19.69423990    50
2   c(0.5, 2, 3, 5, 6) 0.032393618  1.29574473    50
3   c(0.5, 2, 3, 6, 8) 0.473976521 18.95906085    50
4 c(0.5, 1.5, 2, 5, 6) 0.001273863  0.05095453    50

 

 
Associated optimised criterion: 739.2152
 

 
******************* POPULATION FISHER INFORMATION MATRIX ******************
 
                                            s2                        
   173.389525 -43.885540   4.608112 -73.835334   0.000000    0.0000000
   -43.885540  23.487325   1.045124  42.546633   0.000000    0.0000000
     4.608112   1.045124 121.467078   1.590748   0.000000    0.0000000
s2 -73.835334  42.546633   1.590748 148.913217   0.000000    0.0000000
     0.000000   0.000000   0.000000   0.000000 376.204574  295.7085348
     0.000000   0.000000   0.000000   0.000000 295.708535 1036.2908493
     0.000000   0.000000   0.000000   0.000000   2.232828    0.7187292
     0.000000   0.000000   0.000000   0.000000 188.102287  147.8542674
     0.000000   0.000000   0.000000   0.000000 147.854267  518.1454247
     0.000000   0.000000   0.000000   0.000000   1.116414    0.3593646
     0.000000   0.000000   0.000000   0.000000 134.351031   84.0447277
     0.000000   0.000000   0.000000   0.000000 682.199534  -36.5289415
                                                                           
      0.0000000    0.000000 0.000000e+00 0.000000e+00    0.00000    0.00000
      0.0000000    0.000000 0.000000e+00 0.000000e+00    0.00000    0.00000
      0.0000000    0.000000 0.000000e+00 0.000000e+00    0.00000    0.00000
s2    0.0000000    0.000000 0.000000e+00 0.000000e+00    0.00000    0.00000
      2.2328276  188.102287 1.478543e+02 1.116414e+00  134.35103  682.19953
      0.7187292  147.854267 5.181454e+02 3.593646e-01   84.04473  -36.52894
   2952.1468130    1.116414 3.593646e-01 1.476073e+03   85.67843 -586.01956
      1.1164138 3000.629939 4.302075e+03 1.570273e+02  273.35580 1071.43437
      0.3593646 4302.074653 1.020893e+04 6.612112e+01  488.68802 1304.04967
   1476.0734065  157.027322 6.612112e+01 4.687758e+04 2009.31314 5233.27394
     85.6784346  273.355800 4.886880e+02 2.009313e+03 1026.86565 2110.23023
   -586.0195591 1071.434370 1.304050e+03 5.233274e+03 2110.23023 9072.54066

 

 
************************** EXPECTED STANDARD ERRORS ************************
 
------------------------ Fixed Effects Parameters -------------------------
 
                  Beta   StdError       RSE  
ka           1.0000000 0.10506716 10.506716 %
V            3.5000000 0.36476676 10.421907 %
Cl           2.0000000 0.09093069  4.546534 %
beta_V_Sex_2 0.1823216 0.11827387 64.871030 %

 

 
------------------------- Variance of Inter-Subject Random Effects ----------------------
 
   Omega   StdError      RSE  
ka  0.09 0.06744671 74.94079 %
V   0.08 0.03787518 47.34398 %
Cl  0.07 0.01910130 27.28757 %

 
------------------------- Variance of Inter-Occasion Random Effects ----------------------
 
    Gamma    StdError      RSE  
ka 0.0225 0.030090264 133.7345 %
V  0.0125 0.016282630 130.2610 %
Cl 0.0025 0.004924857 196.9943 %

 
------------------------ Standard deviation of residual error ------------------------
 
           Sigma   StdError       RSE  
sig.interA   0.5 0.04522626  9.045252 %
sig.slopeA   0.1 0.01662737 16.627368 %

 
******************************* DETERMINANT ********************************
 
2.662262e+34
 
******************************** CRITERION *********************************
 
739.2152
 

 
***************************** COMPARISON TEST ********************************
 
                  Beta        95 % CI   exp(Beta)       95 % CI
beta_V_Sex_2 0.1823216 [-0.049;0.414]         1.2 [0.952;1.513]

 
Type I error = 0.05
 
             Expected_power Number_subjects_needed (for a given power=0.9)
beta_V_Sex_2       0.338043               176.8715                        

 

 

 
Time difference of 8 secs
sys.self 
    0.22 

 
