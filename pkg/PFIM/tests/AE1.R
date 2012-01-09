require("PFIM")
AE1 <-
    list(PFIMmodel=function(tim, parModel, parArms, model)
     {
         V <- parModel[1]
         k <- parModel[2]
         Alin <- parModel[3]
         dose <- parArms[1]
         nDoses <- parArms[2]
         
         PK <- 0
         for (idx in 0:(nDoses-1)) {
             PK <- PK + (tim>=idx*24)*dose/V*(exp(-k*(tim-idx*24)))
         }
         PD <- Alin*PK
         
         cbind(PK,PD)
     },
         Type='AE',
         parModelName=c('V','k','Alin'),
         parModel=c(4.5,0.5,3),
         parModelVar=c(0.2,0.4,0.15),
         parModelVarType='exp',
         parObsName=c('Conc','Effect'),
         parObsErr=list(c(0.3,0.5), c(0,0.3)),
         parObsTimes=
         list(list(c(0.5,1,2,4,8,12,24,48),
                   c(0,4,8,12,24,48)),
              list(c(0.5,1,2,4,8,12,24),
                   c(0,4,8,12,24)),
              list(c(0.5,1,4,12,23.9,47.9,71.9),
                   c(0,4,12,23.9,47.9,71.9))),
         parArmsName=c('dose','nDoses'),
         parArms=list(c(100,30,5), c(1,1,3)),
         ArmsName=list('100 mg','30 mg','10 mg'),
         TimeName='Time (hr)',
         tRange=c(0,72),
         mpOpt=list()
         )

str(mm <- PFIMmod(AE1))

str(solve(mm))

str(sensitivity(mm))

sensitivity(mm)
