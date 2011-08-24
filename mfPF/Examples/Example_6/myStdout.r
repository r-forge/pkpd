PFIM 3.2 Option 1 
 
Project:  Example 6 
 
Date:  Tue Dec 21 00:05:16 2010
 

 
**************************** INPUT SUMMARY ********************************
 
Analytical function models :  
 
dose/V * ka/(ka - (Cl/V)) * (exp(-(Cl/V) * t) - exp(-ka * t)) 
 
Population design:  
Sample times for response: A 
            times subjects doses
1 c(0.5, 2, 6, 8)       40    30

 
Random effect model: Trand =  2
 
Variance error model response A : ( 0.1 + 0 *f)^2

 
Covariate model :  
 
	NB: Covariates are additive on log parameters
 
		Covariate 1 : Sex ( V ) 
    Categories References Proportions
(1)          M          *         0.5
(2)          F                    0.5

 
Computation of the Fisher information matrix: option =  1
 
******************* POPULATION FISHER INFORMATION MATRIX ******************
 
                                              s2                          
   342.150962 -20.4962991   3.7333850 -30.587703    0.0000000    0.0000000
   -20.496299  31.5727521   0.9836867  57.073989    0.0000000    0.0000000
     3.733385   0.9836867 109.0761038   1.237472    0.0000000    0.0000000
s2 -30.587703  57.0739890   1.2374720 199.758962    0.0000000    0.0000000
     0.000000   0.0000000   0.0000000   0.000000 1463.4668526   65.7219020
     0.000000   0.0000000   0.0000000   0.000000   65.7219020 1871.8838767
     0.000000   0.0000000   0.0000000   0.000000    0.7212424    0.6395265
     0.000000   0.0000000   0.0000000   0.000000  667.6716233  378.8205873
                          
      0.0000000    0.00000
      0.0000000    0.00000
      0.0000000    0.00000
s2    0.0000000    0.00000
      0.7212424  667.67162
      0.6395265  378.82059
   2379.5321132   77.43697
     77.4369654 9002.93885

 
************************** EXPECTED STANDARD ERRORS ************************
 
------------------------ Fixed Effects Parameters -------------------------
 
                  Beta   StdError       RSE  
ka           1.0000000 0.05519603  5.519603 %
V            3.5000000 0.25949173  7.414049 %
Cl           2.0000000 0.09578998  4.789499 %
beta_V_Sex_F 0.4054651 0.10182064 25.112059 %

 
------------------------- Variance of Inter-Subject Random Effects ----------------------
 
   Omega   StdError      RSE  
ka  0.09 0.02660122 29.55691 %
V   0.09 0.02321855 25.79839 %
Cl  0.09 0.02050298 22.78109 %

 
------------------------ Standard deviation of residual error ------------------------ 
 
           Sigma   StdError      RSE  
sig.interA   0.1 0.01076406 10.76406 %

 
******************************* DETERMINANT ********************************
 
6.130894e+21
 
******************************** CRITERION *********************************
 
528.9815
 
***************************** COMPARISON TEST ********************************
 
                  Beta       95 % CI   exp(Beta)       95 % CI
beta_V_Sex_F 0.4054651 [0.206;0.605]         1.5 [1.229;1.831]

 
Type I error = 0.05
 
             Expected_power Number_subjects_needed (for a given power=0.9)
beta_V_Sex_F       0.978421               26.50458                        

 

 

 

 
Time difference of 0.3440001 secs
