PFIM 3.2.1mf  
 
Option: 1 

 
Project: Example 6b 
 
Date: Wed Jun 22 00:28:48 2011
 

 
**************************** INPUT SUMMARY ********************************
 
Analytical function model:  
 
dose/V * ka/(ka - (Cl/V)) * (exp(-(Cl/V) * t) - exp(-ka * t)) 
 

 
Initial population design: 

 
Sample times for response: A 
          Protocol subjects doses
1 c=(0.5, 2, 6, 8)       40    30

 
Total number of samples: 160
 
Associated criterion value: 0
 
Identical sampling times for each response: FALSE
 
Random effect model: Trand =  2
 
Variance error model response A : ( 0.1 + 0 *f)^2

 
Covariate model :  
 
	NB: Covariates are additive on log parameters
 
		Covariate 1 : Sex ( V ) 
    Categories References Proportions
(1)          M          *         0.5
(2)          F                    0.5


Optimization step:  
 
Sampling windows for the response: A 
Window 1 : t= 0.5 2 6 8 
    Nb of sampling points to be taken in this window, n[ 1 ]= 4 
Maximum total number of points in one elementary protocol : 4 
Minimum total number of points in one elementary protocol : 4 

 

Now evaluating the Fisher Information Matrix for the 1 protocols generated 

 
**************************** OPTIMISED DESIGN *****************************
 

 
Optimised population design: 
Sample times for response: A 
            times freq Subjects doses
1 c(0.5, 2, 6, 8)    1       40    30

 

 
Associated optimised criterion: 842.1591
 

 
******************* POPULATION FISHER INFORMATION MATRIX ******************
 
                                         s2                               
   349.67439 -18.46019  211.86402 -28.94599    0.0000    0.0000     0.0000
   -18.46019  32.12770   56.46288  57.45640    0.0000    0.0000     0.0000
   211.86402  56.46288 6050.97836  60.08782    0.0000    0.0000     0.0000
s2 -28.94599  57.45640   60.08782 201.09739    0.0000    0.0000     0.0000
     0.00000   0.00000    0.00000   0.00000 1528.4036   52.7460   650.6747
     0.00000   0.00000    0.00000   0.00000   52.7460 1937.0924   353.9345
     0.00000   0.00000    0.00000   0.00000  650.6747  353.9345 16813.4058

 

 
************************** EXPECTED STANDARD ERRORS ************************
 
------------------------ Fixed Effects Parameters -------------------------
 
                  Beta   StdError        RSE  
ka           1.0000000 0.05515674  5.5156738 %
V            3.5000000 0.25772605  7.3636014 %
Cl           2.0000000 0.01318052  0.6590262 %
beta_V_Sex_2 0.4054651 0.10103458 24.9181929 %

 

 
------------------------- Variance of Inter-Subject Random Effects ----------------------
 
   Omega   StdError      RSE  
ka  0.09 0.02579898 28.66553 %
V   0.09 0.02277067 25.30075 %

 
------------------------ Standard deviation of residual error ------------------------
 
           Sigma    StdError      RSE  
sig.interA   0.1 0.007789801 7.789801 %

 
******************************* DETERMINANT ********************************
 
3.004408e+20
 
******************************** CRITERION *********************************
 
842.1591
 

 
***************************** COMPARISON TEST ********************************
 
                  Beta       95 % CI   exp(Beta)       95 % CI
beta_V_Sex_2 0.4054651 [0.207;0.603]         1.5 [1.231;1.828]

 
Type I error = 0.05
 
             Expected_power Number_subjects_needed (for a given power=0.9)
beta_V_Sex_2       0.979972               26.09692                        

 

 

 
Time difference of 0.2809999 secs
sys.self 
       0 

 
