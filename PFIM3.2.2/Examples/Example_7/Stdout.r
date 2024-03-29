PFIM 3.2 Option 1 
 
Project:  Example 7 
 
Date:  Thu Feb 10 14:56:30 2011
 

 
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
	Covariate 2 : Genetics ( V ) 
    Categories References Proportions
(1)  common_Hz          *        0.50
(2)         hz                   0.25
(3)    rare_hz                   0.25

 
Computation of the Fisher information matrix: option =  1
 
******************* POPULATION FISHER INFORMATION MATRIX ******************
 
                                             s2          s2          s2                                                 
   341.583768 -21.126433   3.8060206 -37.902510 -18.9016959 -19.1722946    0.0000000    0.000000    0.0000000    0.00000
   -21.126433  31.054082   1.0187123  53.585645  26.7276369  26.7093331    0.0000000    0.000000    0.0000000    0.00000
     3.806021   1.018712 108.7748395   1.802899   0.8870889   0.9569525    0.0000000    0.000000    0.0000000    0.00000
s2 -37.902510  53.585645   1.8028994 187.549757  43.9992742  43.6710020    0.0000000    0.000000    0.0000000    0.00000
s2 -18.901696  26.727637   0.8870889  43.999274  93.5467291   0.0000000    0.0000000    0.000000    0.0000000    0.00000
s2 -19.172295  26.709333   0.9569525  43.671002   0.0000000  93.4826660    0.0000000    0.000000    0.0000000    0.00000
     0.000000   0.000000   0.0000000   0.000000   0.0000000   0.0000000 1459.0294554   70.387310    0.7592298  661.52089
     0.000000   0.000000   0.0000000   0.000000   0.0000000   0.0000000   70.3873104 1813.866409    0.6891870  411.22779
     0.000000   0.000000   0.0000000   0.000000   0.0000000   0.0000000    0.7592298    0.689187 2366.4460382   88.78356
     0.000000   0.000000   0.0000000   0.000000   0.0000000   0.0000000  661.5208877  411.227789   88.7835606 9081.52067

 
************************** EXPECTED STANDARD ERRORS ************************
 
------------------------ Fixed Effects Parameters -------------------------
 
                             Beta   StdError       RSE  
ka                      1.0000000 0.05530278  5.530278 %
V                       3.5000000 0.31461129  8.988894 %
Cl                      2.0000000 0.09592383  4.796191 %
beta_V_Sex_F            0.4054651 0.10267337 25.322370 %
beta_V_Genetics_hz      0.2623643 0.12603290 48.037372 %
beta_V_Genetics_rare_hz 0.3364722 0.12607617 37.470006 %

 
------------------------- Variance of Inter-Subject Random Effects ----------------------
 
   Omega   StdError      RSE  
ka  0.09 0.02663195 29.59106 %
V   0.09 0.02360907 26.23230 %
Cl  0.09 0.02056052 22.84502 %

 
------------------------ Standard deviation of residual error ------------------------ 
 
           Sigma   StdError      RSE  
sig.interA   0.1 0.01072185 10.72185 %

 
******************************* DETERMINANT ********************************
 
2.532051e+25
 
******************************** CRITERION *********************************
 
347.0142
 

 

 

 

 
Time difference of 0.06599998 secs
