require("PFIM")
AEfun3 <-
    list(PFIMmodel=function(tim, parModel, parArms, model)
     {
         V <- parModel[1]
         k <- parModel[2]
         keff <- parModel[3]
         Emax <- parModel[4]
         EC50 <- parModel[5]
         dose <- parArms[1]
         
         if (tim <= 7) { # Appearance of ADAs at t=7
             PK = dose/V*(exp(-k*tim))
             temp = dose/V*k/(k-keff)*(exp(-keff*tim)-exp(-k*tim))
         } else {
             PK = 0
             temp = 0
         }
         Eff <- Emax*temp/(EC50 + temp)
         PD <- 1+Eff
         
         return(cbind(PK,PD))
     },
         Type='AE',
         parModelName=c('V','k','keff','Emax','EC50'),
         parModel=c(4.5,0.15,3,2.3,3),
         parModelVar=c(0.2,0.4,0.15,0.4,0.05),
         parModelVarType='exp',
         parObsName=c('Conc','Effect'),
         parObsErr=list(c(0.3,0.5), c(0.2,0.7)),
         ## Sampling times per arm - and within arm per observation/measurement
         parObsTimes=list(c(0.5,1,2,4,8,12,24,48)),
         parArmsName=c('dose'),
         parArms =list(c(200)),
         ArmsName=list('PKPD model'),
         TimeName='Time (hr)',
         tRange=c(0,48),
         mpOpt=list()
     )

str(mm <- PFIMmod(AEfun3))
str(solve(mm))
