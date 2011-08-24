PFIM 3.2  
 
Option: 1 

 
Project: Example 9
 
Date: Tue Dec 21 00:26:02 2010
 

 
**************************** INPUT SUMMARY ********************************
 
Analytical function model:  
 
dose/V * ka/(ka - (Cl/V)) * (exp(-(Cl/V) * t) - exp(-ka * t)) 
 

 
Initial population design: 

 
Sample times for response: A 
             Protocol subjects doses
1 c=(0.5, 2, 4, 6, 8)       20    30
2   c=(1, 2, 3, 4, 5)       20    50

 
Total number of samples: 200
 
Associated criterion value: 1590.883
 
Identical sampling times for each response: FALSE
 
Number of occasions: 2
 
Random effect model: Trand =  2
 
Variance error model response A : ( 0.1 + 0 *f)^2

 
Covariate model :  
 
	NB: Covariates are additive on log parameters
 
	Covariates not changing with occasion 
 
	Covariate 1 : Sex ( V ) 
    Categories References Proportions
(1)          M          *         0.5
(2)          F                    0.5


 
	Covariates changing with occasion 
 
	Covariate  1 : Treat ( Cl ) 
    Categories References
(1)          A          *
(2)          B           

 
    Sequences Proportions
(1)       A B         0.5
(2)       B A         0.5

 

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
               times freq Subjects doses
1 c(0.5, 2, 5, 6, 8)    1       40    50

 

 
Associated optimised criterion: 1740.182
 

 
******************* POPULATION FISHER INFORMATION MATRIX ******************
 
                                              s2                         
   372.9169856 -4.7752523  0.5901490  -7.1926884   0.8510045 0.000000e+00
    -4.7752523 31.1486508  0.1340608  54.8441354   0.2238203 0.000000e+00
     0.5901490  0.1340608 98.5844358   0.2090275  98.4875256 0.000000e+00
s2  -7.1926884 54.8441354  0.2090275 191.9544739   0.3079808 0.000000e+00
     0.8510045  0.2238203 98.4875256   0.3079808 973.3587964 0.000000e+00
     0.0000000  0.0000000  0.0000000   0.0000000   0.0000000 1.738388e+03
     0.0000000  0.0000000  0.0000000   0.0000000   0.0000000 3.559460e+00
     0.0000000  0.0000000  0.0000000   0.0000000   0.0000000 1.745120e-02
     0.0000000  0.0000000  0.0000000   0.0000000   0.0000000 8.707543e+02
     0.0000000  0.0000000  0.0000000   0.0000000   0.0000000 3.332722e+00
     0.0000000  0.0000000  0.0000000   0.0000000   0.0000000 1.081155e-02
     0.0000000  0.0000000  0.0000000   0.0000000   0.0000000 4.005766e+02
                                                                               
   0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00     0.00000
   0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00     0.00000
   0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00     0.00000
s2 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00     0.00000
   0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00     0.00000
   3.559460e+00 1.745120e-02 8.707543e+02 3.332722e+00 1.081155e-02   400.57662
   1.820023e+03 1.113889e-02 3.054206e+00 9.112824e+02 8.374689e-03   243.57056
   1.113889e-02 1.943778e+03 1.146938e-02 9.767236e-03 9.718900e+02    14.33021
   3.054206e+00 1.146938e-02 2.576895e+04 2.023354e+03 9.828009e+00  2930.71816
   9.112824e+02 9.767236e-03 2.023354e+03 3.378150e+04 6.483263e+00  1829.03850
   8.374689e-03 9.718900e+02 9.828009e+00 6.483263e+00 4.739027e+04   270.57813
   2.435706e+02 1.433021e+01 2.930718e+03 1.829039e+03 2.705781e+02 33847.86827

 

 
************************** EXPECTED STANDARD ERRORS ************************
 
------------------------ Fixed Effects Parameters -------------------------
 
                      Beta   StdError       RSE  
ka              1.00000000 0.05183605  5.183605 %
V               3.50000000 0.25433707  7.266774 %
Cl              2.00000000 0.10622799  5.311400 %
beta_V_Sex_2    0.18232156 0.10239047 56.159279 %
beta_Cl_Treat_B 0.09531018 0.03380680 35.470289 %

 

 
------------------------- Variance of Inter-Subject Random Effects ----------------------
 
   Omega   StdError      RSE  
ka  0.09 0.02421047 26.90052 %
V   0.09 0.02360872 26.23191 %
Cl  0.09 0.02279896 25.33218 %

 
------------------------- Variance of Inter-Occasion Random Effects ----------------------
 
    Gamma    StdError      RSE  
ka 0.0225 0.006323905 28.10625 %
V  0.0225 0.005496691 24.42974 %
Cl 0.0225 0.004617460 20.52204 %

 
------------------------ Standard deviation of residual error ------------------------
 
           Sigma    StdError      RSE  
sig.interA   0.1 0.005474988 5.474988 %

 
******************************* DETERMINANT ********************************
 
7.711473e+38
 
******************************** CRITERION *********************************
 
1740.182
 

 

 

 
***************************** EQUIVALENCE TEST ********************************
 
                      Beta       90 % CI   exp(Beta)       90 % CI
beta_V_Sex_2    0.18232156 [0.014;0.351]         1.2  [1.014;1.42]
beta_Cl_Treat_B 0.09531018  [0.04;0.151]         1.1 [1.041;1.163]

 
Type I error = 0.05
Equivalence interval = [log(0.8),log(1.25)]
 
                Expected_power Number_subjects_needed (for a given power=0.9)
beta_V_Sex_2         0.1063521             2155.06096                        
beta_Cl_Treat_B      0.9836782               23.95788                        

 
Time difference of 17.562 secs
sys.self 
    0.25 

 
