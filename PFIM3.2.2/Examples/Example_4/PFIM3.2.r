#PFIM 3.2.2 Multiple responses
#February 2011
#Copyright © PFIM 3.2.2 – Caroline Bazzoli, Thu Thuy Nguyen, Anne Dubois, Sylvie Retout, Emanuelle Comets, France Mentré - Université Paris Diderot – INSERM.
#------------------------------------------------------------------------------------

rm(list=ls(all=TRUE))

under.unix<-!(version$os=='Microsoft Windows' ||
      version$os=='Win32' || version$os=='mingw32')
dirsep<-ifelse(under.unix,"/","\\")

library(nlme)
library(odesolve)
options(expressions=10000,max.deparse.lines=10000,max.deparse.length=10000,deparse.max.length=10000)

subjects<-NULL

directory<-"C:\\Users\\Thu Thuy Nguyen\\Desktop\\PFIM3.2.2\\Examples\\Example_4"
directory.program<-"C:\\Users\\Thu Thuy Nguyen\\Desktop\\PFIM3.2.2\\PFIM3.2.2\\Program"


#############################

PFIM<-function(model.file="stdin.r") {
   under.unix<-!(version$os=='Microsoft Windows' ||
      version$os=='Win32' || version$os=='mingw32')
   dirsep<-ifelse(under.unix,"/","\\")
   cat("Directory for PFIM3.2 programs:",directory.program,"\n")
   cat("Directory for project:",directory,"\n")
   source(paste(directory,dirsep,model.file,sep=""))
   source(paste(directory,dirsep,file.model,sep=""))
     # download the file "modelcreated.r" linked to the R function create_formED()
    
    if(length(grep("CreateModel_PKPDdesign.r",readLines(paste(directory,dirsep,file.model,sep=""),warn = F)))!=0){
     source(paste(directory,dirsep,"model_created.r",sep=""))}

  #test to accept stdin file of version PFIM 3.0

    #names.datax and names.datay objects
        ngraph<-2
        vec<-c("x","y")
        err1<-tryCatch(names.data_test<-lapply(1:ngraph,function(ngraph,x,vec) get(x), x=paste("names.data",vec[ngraph],sep=""), vec=vec), error=function(e) 4)

         if(length(which(is.null(unlist(err1))))>=1 || length(which(is.null(unlist(err1)==4)))>=1)
         {
         names.datax<-rep("Time",nr)
         names.datay<-rep("Concentration",nr)
         }


    #option object
      err2<-tryCatch(get("option"), error=function(e) 4)
      if(err2==4) option<-1

    #covariate.model
      err3<-tryCatch(get("covariate.model"), error=function(e) 4)
      if(err3==4) covariate.model<-F

    #IOV
      err4<-tryCatch(get("n_occ"), error=function(e) -4)
      if(err4==-4) n_occ<-0

     #covariate_occ.model
        err5<-tryCatch(get("covariate_occ.model"), error=function(e) 4)
        if(err5==4) covariate_occ.model<-F



     if (run=="EVAL"){
        source(paste(directory,file.model,sep=dirsep))
        if (modelform=="DE")
        {
          if(option==2)
          {
            if (covariate.model==T || n_occ>1) {stop("You can not use option 2 with covariate or inter-occasion variability currently")}
            else {source(paste(directory.program,dirsep,"EQPfim3.2op2.r",sep=""))}
            out()
           }
           else
           {

            source(paste(directory.program,dirsep,"EQPfim3.2op1.r",sep=""))
            out()

           }

        }
        else
        {
          if(option==2)
          {
            if (covariate.model==T || n_occ>1) {stop("You can not use option 2 with covariate or inter-occasion variability currently")}
            else {source(paste(directory.program,dirsep,"Pfim3.2op2.r",sep="")) }
            out()
          }
          else
          {

             source(paste(directory.program,dirsep,"Pfim3.2op1.r",sep=""))
             out()

           }
          }
        }

          else {

            if (algo.option=="FW")
                {    #if (n_occ>1) {stop("You can not use Fedorov-Wynn algorithm with inter-occasion variability currently")}
                     source(paste(directory.program,dirsep,"algofedorov3.2.r",sep=""))
#provisoire : pour l'instant je source rcostfed.pr, mais il devrait être renommé
# version normale de fedorov rdosefed.pr (à renommer en algofedorov.r par ex)
  #          source(paste(directory.program,"\\algofedorov.simu.r",sep=""))
            out.fedorov(modelfile=model.file,directory=directory,
               directory.program=directory.program)
                }

            else
            {


                if (modelform=="DE")
                {
                   if (option==1)
                    source(paste(directory.program,dirsep,"EQPfimOPT3.2op1.r",sep=""))
                    else {
                    if (covariate.model==T || n_occ>1) {stop("You can not use option 2 for Simplex algorithm with covariate or inter-occasion variability currently")}
                    else {source(paste(directory.program,dirsep,"EQPfimOPT3.2op2.r",sep=""))}
                    }
                }
                else
                {
                    if (option==1)
                    source(paste(directory.program,dirsep,"PfimOPT3.2op1.r",sep=""))
                     else {
                    if (covariate.model==T || n_occ>1) {stop("You can not use option 2 for Simplex algorithm with covariate or inter-occasion variability currently")}
                    else {source(paste(directory.program,dirsep,"PfimOPT3.2op2.r",sep=""))}
                    }
                }

              if (algo.option=="SIMP")
              {
               source(paste(directory,file.model,sep=dirsep))
               source(paste(directory.program,dirsep,"algosimplex3.2.r",sep=""))
               x<-out.simplex()
		           print(x)
               #remove(subjects,inherits = FALSE)
              }
              else
	             cat("Error in the specification of the optimisation algorithm.","\n")
          }
     }
}