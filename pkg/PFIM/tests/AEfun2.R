require("PFIM")
AEfun2 <-
     list(PFIMmodel=function(tim, parModel, parArms, model)
      {
          V <- parModel[1]
          k <- parModel[2]
          keff <- parModel[3]
          Emax <- parModel[4]
          EC50 <- parModel[5]
          dose <- parArms[1]
          nDoses <- parArms[2]

          PK <- 0
          temp <- 0
          II <- 12
          for (idx in 0:(nDoses-1)) {
              PK <- PK + (tim>=idx*II)*dose/V*(exp(-k*(tim-idx*II)))
              temp <- temp + (tim>=idx*II)*dose/V*k/(k-keff)*(exp(-keff*(tim-idx*II))-exp(-k*(tim-idx*II)))
          }
          Eff <- Emax*temp/(EC50 + temp)
          PD <- exp(-Eff)

          return(cbind(PK,temp,PD))
      },
          Type='AE',
         parModelName=c('V','k','keff','Emax','EC50'),
         parModel=c(4.5,0.15,3,120,3),
         parModelVar=c(0.2,0.4,0.15,0.4,0.05),
         parModelVarType='exp',
         parObsName=c('Conc','Effect 1','Effect 2'),
         parObsErr=list(c(0.3,0.5), c(0,0.3), c(0.2,0.7)),
         # Sampling times per arm - and within arm per observation/measurement
         parObsTimes=
          list(list(
                    c(0.5,1,2,4,8,12,24,48),
                    c(0,4,8,12,24,48),
                    c(-0.1,0,4,8,12,24,48)),
               list(c(0.5,1,2,4,8,11.9,23.9),
                    c(0,4,8,12,23.9,47.9),
                    c(-0.1,0,4,8,12,23.9,47.9))),
         parArmsName=c('dose','nDoses'),
         parArms =list(c(30,5), c(1,3)),
         ArmsName=list('30 mg','10 mg'),
         TimeName='Time (hr)',
         tRange=c(-1,48),
         mpOpt=list()
      )
str(mm <- PFIMmod(AEfun2))

str(solve(mm))
