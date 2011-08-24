#PFIM 3.2.2 Multiple responses
#February 2011
#Copyright © PFIM 3.2.2 – Caroline Bazzoli, Thu Thuy Nguyen, Anne Dubois, Sylvie Retout, Emanuelle Comets, France Mentré - Université Paris Diderot – INSERM.
#------------------------------------------------------------------------------------


                                                                     

#Version R : modification par rapport à Splus
#1- ne pas mettre de '; '
#2- retirer les 'break '                                                                       
#3- remplacer 'null' par 'NULL'
#4- cat(paste(form[..]...) à la place de cat(form[...]...)
#5- ne pas couper la phrase dans le cat('Error:...")et laisser le ';'
#6- dans le Simplex, les "cat" suivis des "print" doivent se trouver sur deux lignes differentes
#7- changer la fonction determinant en det.
###############


###########################################################################
# 	Optimisation of population designs from computation of the incomplete 
#	expression of the Population Fisher Information matrix for PK Design.     	           
#	Use of a Simplex algorithm.							 					     
###########################################################################


  

#---------------------- PRE - COMPUTATION ------------------------------------
#-----------------------------------------------------------------------------

#test to accept stdin file of version PFIM 3.0
#option object
err2<-tryCatch(get("option"), error=function(e) 4)
if(err2==4) option<-1

#names.datax and names.datay objects
ngraph<-2
vec<-c("x","y")
err1<-tryCatch(names.data_test<-lapply(1:ngraph,function(ngraph,x,vec) get(x), x=paste("names.data",vec[ngraph],sep=""), vec=vec), error=function(e) 4)
 if(length(err1)==1 && err1==4)
 {
 names.datax<-rep("Time",nr)
 names.datay<-rep("Concentration",nr)
 }
 
#covariate.model
err3<-tryCatch(get("covariate.model"), error=function(e) 4)
if(err3==4) covariate.model<-F

#IOV
err4<-tryCatch(get("n_occ"), error=function(e) -4)
if(err4==-4) n_occ<-0

#covariate_occ.model
err5<-tryCatch(get("covariate_occ.model"), error=function(e) 4)
if(err5==4) covariate_occ.model<-F
    
  
 
#----------------------------------------

  theta<-beta
  lf<-length(form)#number of analytical form
#nombre de form dans chaque réponse pour savoir quoi prendre en tf
#recover of the model
#partie liée à l'implémentation pour reponses multiples
  formg<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("form",LETTERS[1:nr],sep=""))  #sous forme de list
  formg_init<-formg
  lformg<-lapply(formg,length)  #length(of each model tf or not)
  c<-c(0,cumsum(lformg))
  doses<-dose
  p<-length(theta)
#Partie liée à l'implémentation pour iov
#-------------------------------------------------------------------------------
if (n_occ<=1) 
{
  n_occ=1
  covariate_occ.model<-F
}

if (covariate_occ.model==F)
{
  covariate_occ.name<-list(NULL)
  covariate_occ.category<-list(NULL)
  covariate_occ.sequence<-list(list(NULL))
  covariate_occ.proportions<-list(NULL)
  parameter_occ.associated<-list(NULL)
  beta.covariate_occ<-list(NULL)
  locc<-0
  nb_seq<-1
}
sequence.proportions<-covariate_occ.proportions

 
pn<-n_occ
v<-rep(0,p)
param.asso.occ<-list()
occ.effect<-list()
beta.occ<-list()
nb_seq<-c()
  
for (kj in 1:length(covariate_occ.name)){
    nb_seq[kj]=length(covariate_occ.sequence[[kj]])
    param.asso.occ[[kj]]<-list()
    occ.effect[[kj]]<-list()
    beta.occ[[kj]]<-list()
    for (l in 1:nb_seq[kj]) {
      param.asso.occ[[kj]][[l]]<-list()
      occ.effect[[kj]][[l]]<-list()
      beta.occ[[kj]][[l]]<-list()
      for (i in 1:n_occ)  {param.asso.occ[[kj]][[l]][[i]]=v;occ.effect[[kj]][[l]][[i]]=v;beta.occ[[kj]][[l]][[i]]=v} 
      }
     
    if (n_occ>1){
      for (l in 1:nb_seq[kj]) {
  
        num<-c()
        if (covariate_occ.model==T){
          for (i in which(covariate_occ.sequence[[kj]][[l]]!=covariate_occ.category[[kj]][1])){
            num[kj]=1
            param.asso.occ[[kj]][[l]][[i]]=parameter_occ.associated[[kj]]
            
            for (j  in 1:length(param.asso.occ[[kj]][[l]][[i]])){
              
              occ.effect[[kj]][[l]][[i]][which(param.asso.occ[[kj]][[l]][[i]][j]==parameters)]=1
              beta.occ[[kj]][[l]][[i]][which(param.asso.occ[[kj]][[l]][[i]][j]==parameters)]=beta.covariate_occ[[kj]][[j]][num[kj]]
   
            } 
            num[kj]=num[kj]+1
         }
       }
     }
  }
}

#les combinaisons des séquences
sequence.reco<-list()
for (kj in 1:length(covariate_occ.name)){
  sequence.reco[[kj]]=c(1:nb_seq[kj])
}
comb.seq.poss<-unique(t(combn(unlist(sequence.reco),length(covariate_occ.name))))
ll.cov_occ<-length(covariate_occ.name)
for (i in 1: ll.cov_occ)
{
  nline<-which(comb.seq.poss[,i] > max(sequence.reco[[i]]) | comb.seq.poss[,i] < min(sequence.reco[[i]]))
  ref<-integer(length = 0)
  if (!identical(nline,ref)) comb.seq.poss<-comb.seq.poss[-nline,]
}
#nb de sequences possibles          
lposs<-dim(comb.seq.poss)[1]

#les propportions de chaque combinaison de séquence
comb.seq.prop<-c()
nbcov_occ<-length(covariate_occ.name)
if (covariate_occ.model==T){
for (i in 1: lposs){
comb.seq.prop[i]<-sequence.proportions[[1]][comb.seq.poss[i,1]]
if (nbcov_occ>1){
  for (j in 2:nbcov_occ){
    comb.seq.prop[i]=comb.seq.prop[i]*sequence.proportions[[j]][comb.seq.poss[i,j]]
  }
 }
}
}else{comb.seq.prop<-1}

#le vecteur des parametres fixes prenant compte des effets covariables changeant entre les périodes

theta1<-list()
Pt1<-list()
vecP1<-list()
lP1<-list() 
for (kj in 1:length(covariate_occ.name)){
  theta1[[kj]]<-list()
  Pt1[[kj]]<-list()
  vecP1[[kj]]<-list()
  lP1[[kj]]<-list()
  for (l in 1:nb_seq[kj]) {
  theta1[[kj]][[l]]<-list()
  Pt1[[kj]][[l]]<-list()
  vecP1[[kj]][[l]]<-list()
  lP1[[kj]][[l]]<-rep(0,pn)
  for (i in 1:pn)  { 
    Pt1[[kj]][[l]][[i]]<-c(occ.effect[[kj]][[l]][[i]])
    vecP1[[kj]][[l]][[i]]<-Pt1[[kj]][[l]][[i]]!=0
    Pt1[[kj]][[l]][[i]]<-Pt1[[kj]][[l]][[i]][vecP1[[kj]][[l]][[i]]]
    lP1[[kj]][[l]][i]<-length(Pt1[[kj]][[l]][[i]])
    if (lP1[[kj]][[l]][i] == 0) {
    theta1[[kj]][[l]][[i]] = beta
    }
    if (lP1[[kj]][[l]][i] > 0){
      if (Trand==1) theta1[[kj]][[l]][[i]]<-beta+beta.occ[[kj]][[l]][[i]]
      if (Trand==2) theta1[[kj]][[l]][[i]]<-beta*exp(beta.occ[[kj]][[l]][[i]])
    }
  }
 }
}

vecP<-list()
beta.occ2<-list()
lP<-list() 
for (i in 1: lposs){
 vecP[[i]]<-list()
 lP[[i]]<-rep(0,pn)
 beta.occ2[[i]]<-list()
 for (occ in 1: n_occ){
  vecP[[i]][[occ]]<-vecP1[[1]][[comb.seq.poss[i,1]]][[occ]]
  beta.occ2[[i]][[occ]]<-beta.occ[[1]][[comb.seq.poss[i,1]]][[occ]]
  lP[[i]]<-lP1[[1]][[comb.seq.poss[i,1]]]
  if (nbcov_occ>1){
   for (j in 2:nbcov_occ){
    vecP[[i]][[occ]]=vecP[[i]][[occ]]*vecP1[[j]][[comb.seq.poss[i,j]]][[occ]]
    beta.occ2[[i]][[occ]]=beta.occ2[[i]][[occ]]+beta.occ[[j]][[comb.seq.poss[i,j]]][[occ]]
    lP[[i]]<-lP[[i]]+lP1[[j]][[comb.seq.poss[i,j]]]
   }
  }
 }
}


theta2<-list()
  for (l in 1:lposs) {
  theta2[[l]]<-list()
  for (i in 1:n_occ)  { 
    if (lP[[l]][i] == 0) {
    theta2[[l]][[i]] = beta
    }
    if (lP[[l]][i] > 0){
      if (Trand==1) theta2[[l]][[i]]<-beta+beta.occ2[[l]][[i]]
      if (Trand==2) theta2[[l]][[i]]<-beta*exp(beta.occ2[[l]][[i]])
    }
  }
 }

#valeurs des effets covariables qui dependent de l'occasion  
beta.occ1<-unlist(beta.covariate_occ)    
locc<-length(beta.occ1)

#pour prendre en compte les variances des effets aléatoires liés à IOV

if (pn==1) gamma<-diag(NULL)
omega<-as.matrix(omega)
gamma<-as.matrix(gamma)
OMEGA<-diag(c(diag(omega),rep(c(diag(gamma)),pn))) #la matrice OMEGA qui regroupe omega et gamma
if (length(c(diag(omega),rep(c(diag(gamma)),pn)))==1) OMEGA<-as.matrix((c(diag(omega),rep(c(diag(gamma)),pn))))

#----------- Fin de la partie precomputation pour IOV -------------------------------------------



#recover the inital data for all the models
  prots<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("prot",LETTERS[1:nr],sep=""))
  sigmainter<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("sig.inter",LETTERS[1:nr],sep=""))
  sigmaslope<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("sig.slope",LETTERS[1:nr],sep=""))

  #tests if all bound are written
  pres_b<-lapply(1:nr,function(nr,x) exists(x[nr]),x=paste("bound",LETTERS[1:nr],sep="")) 
  TF<-which(pres_b==F) 
  if(length(TF)==0) bornes<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("bound",LETTERS[1:nr],sep="")) else {

  if(length(TF)==nr)  bornes<-rep(list(list(c(0,Inf))),nr) else { 
    bornes<-rep(list(list(c(0,Inf))),nr)
    for (i in 1:nr){ 
          if(length(which(TF==i))==0 ) {
             x=paste("bound",LETTERS[i],sep="")
             bornes[[i]]<-get(x)
          }
    }
 }
}
       
#test if graph.inf existe pour chaque réponse#
erro<-tryCatch(graphinf<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("graph.inf",LETTERS[1:nr],sep="")),error=function(e) 1)          
#e=0 si pas d'erreur et 1 sinon#
if (length(which(unlist(erro)==1))>=1)
{
vecg<-replicate(nr,0)
graphinf<- as.vector(vecg,mode="list")
graphsup<-lapply(1:nr,function(nr,x) max(unlist(get(x[nr]))),x=paste("prot",LETTERS[1:nr],sep=""))
} else
{
graphinf<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("graph.inf",LETTERS[1:nr],sep=""))
graphsup<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("graph.sup",LETTERS[1:nr],sep=""))
}

#test si y.range exist y.rangeA<-NULL # default range
erro1<-tryCatch(y.range<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("y.range",LETTERS[1:nr],sep="")),error=function(e) 1)          
if (length(which(unlist(erro1)==1))>=1)
{
vecg1<-replicate(nr,NULL)
y.range<- as.vector(vecg1,mode="list")
} else
{
y.range<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("y.range",LETTERS[1:nr],sep=""))
}


#recover tf in order to compute the sensibility
tf<-vector(mode="list")
for (i in 1 : length(lformg)){

  if (lformg[i]!=1){ 
    tf[[i]]<-vector()
    tf[[i]][1]<-bornes[[i]][[1]][2]
    for (j in 2:lformg[[i]]){
      
      tf[[i]][j]<-c(bornes[[i]][[j]][2]) 
    }
  }
  else {
      tf[[i]]<-bornes[[i]][[1]][2]
       #tf[[i]]<-bornes[[i]][2]
  }     
} 

#allow to build graph 
tt1<-list()
		#tt<-lapply(1:nr,function(nr,graphinf,graphsup) seq(graphinf,graphsup,0.1),graphinf<-graphinf[[nr]],graphsup<-graphsup[[nr]])	
	  for (l in 1:nr){
          tt1[[l]]<-list()
          if (length(bornes[[l]])==1){
            tt1[[l]][[1]]<-c(bornes[[l]][[1]][1]:graphsup[[l]])
            
           } else {
  		      for (i in 1: (length(bornes[[l]])-1)){
  		        tt1[[l]][[i]]<-c(bornes[[l]][[i]][1]:bornes[[l]][[i]][2])
  		        if(i==length(bornes[[l]])-1){
  		         # if (bornes[[l]][[i]][2]>graphsup[[l]]) stop(" The upper limit (graph.sup)  for the response  ",LETTERS[ l]  ,"  is too small")
  		          tt1[[l]][[i+1]]<-c(bornes[[l]][[i]][2]:graphsup[[l]])}
             }
          }
    }
tt1b<-lapply(1:nr,function(nr,tt1) unlist(tt1[[nr]]),tt1<-tt1) 
#tt1<-lapply(1:nr,function(nr,tt1b) tt1b[[nr]][tt1b[[nr]]>=graphinf[[nr]][1] & tt1b[[1]]<=graphsup[[nr]][1]],tt1b<-tt1b) 
tt1<-lapply(1:nr,function(nr,tt1b) tt1b[[nr]][tt1b[[nr]]>=graphinf[[nr]][1] & tt1b[[nr]]<=graphsup[[nr]][1]],tt1b<-tt1b) 

#vector of dose
if (dose.identical==T) doses<-rep(doses,length(prots[[1]])) else NULL

#evaluation of individual design 
b_condition<-T
if (length(omega)!=0 && length(which(omega==0))==length(omega)) 
{b_condition<-F;
Omega<-NULL;
cat("All elements of the vector omega are equal to 0, no random effect is considered","\n")
#ls<-length(subjects); subjects<-rep(1,ls)
}



#-----------------------------------------------------------------------------
#------------------------ COMPUTE --------------------------------------------
#Function used in Population Fisher information matrix  function 
#------------------------------------
pts<-function(nr,i,prots){
protk<-c()
lprotk<-c()
for (kj in 1:nr){
protk <-c(protk,prots[[kj]][i])  
}
return(protk)
}

#--------------------------------------------------------------------------------
#To associate proportions of categories in a same tab
#--------------------------------------------------------------------------------
association_modality_covariate<-function(covariate.name,covariate.cat.reco,covariate.proportions,covariate.cat.poss)
{
pro.cov<-NULL
summary.pro.cov<-NULL

lcov.name<-length(covariate.name) 
lmod.poss<-dim(covariate.cat.poss)[1]
prop.comb<-matrix(nrow=lmod.poss,ncol=lcov.name)
for ( i in 1:lcov.name)
{
    pro.cov<-cbind(rep(covariate.name[[i]],length(covariate.cat.reco[[i]])),as.numeric(covariate.cat.reco[[i]]),as.numeric(covariate.proportions[[i]]) )
    #pro.cov<-cbind(covariate.category[[i]],covariate.proportions[[i]])
    summary.pro.cov<-rbind(summary.pro.cov,pro.cov)
}

for (i in 1:lcov.name){
                                          for (ij in 1: lmod.poss)
                                          {
                                          names.cov.ass<-covariate.name[[i]]
                                          mod.cov.ass<-covariate.cat.poss[ij,i]
                                          stab<-summary.pro.cov[summary.pro.cov[,1]==names.cov.ass,]
                                          prop.comb[ij,i]<-as.numeric(stab[stab[,2]==mod.cov.ass,3])
                                          }
                 }                                       

return(list(summary.pro.cov,prop.comb))
}


#Sensitivity matrix arguments:
#nr: number of responses
#protk: elementary design
#pfixe: number of fixed effects
#p: length of the vector beta in the stdin.r file
#kjk: number of elementray design
#lprotk: number of samples in the kith elementray design
#doses
#arguments if covariates
#categories.cov: categories for each covariate
#names.cov: name of each covariate
#beta.covaraite:
#theta: value of fixed effects in the stind.r file
#tab.cov.param
#--------------------
sensibilite<-function(nr,protk,pfixe,p,kjk,lprotk,doses,categories.cov,names.cov,beta.covariate,thetah,tab.cov.param,beta.occh)
{	
	#c<-length(proti)
  dose<-doses[kjk] 
  theta_init<-thetah
  s1<-matrix(nrow=p,ncol=0)
  s4<-matrix(nrow=p,ncol=0)
  #s1rand<-matrix(nrow=prand,ncol=0)
	mod<-numeric()
	mod1<-numeric()
	
	#IF COVARIATE  
	if (!is.null(categories.cov))
	{
        	parameters.cov.val<-rep(1,length(thetah))
        	n.theta<-c(theta_init,as.numeric(unlist(tab.cov.param[,5]))) 
          sum<-0
          for (i in 1:pfixe) {assign(c(parameters,tab.cov.param[,3])[i],n.theta[i])}
        	#attibuer les valeurs des modalités aux covariables
        
        	for (i in 1:length(names.cov)) {assign(names.cov[[i]],categories.cov[i])} #assign("trt",groupe)
     
          #boucle sur le nombre de parametres ayant une cov associée
          for (i in 1: length(unique(tab.cov.param[,1])))
        	{
          	   
            	 par.test<-unique(tab.cov.param[,1])[i]
            	 tab.cov.param1<-tab.cov.param[which(tab.cov.param[,1]==par.test),]
            	 tab.cov.param1<-matrix(tab.cov.param1,ncol=5)
             
             if (length(unique(tab.cov.param1[,2]))==1)
             { 
              # cas avec d'un parametre associé à une covariable  à une ou plusiers modalités
              # ex (Trand==2): Cl*exp(beta_CL*COV)
              #(Trand==1):Cl + beta_CL*COV
                test.bet<-tab.cov.param1[get(tab.cov.param1[1,2])==tab.cov.param1[,4],3]
                
                if (length(test.bet)!=0)  #si categories ==0
                {
                  if (Trand==2)   
                  {
                     
                      thetah[which(tab.cov.param1[1,1]==parameters)]<-thetah[which(tab.cov.param1[1,1]==parameters)]*exp(get(tab.cov.param1[get(tab.cov.param1[1,2])==tab.cov.param1[,4],3]))
                      parameters.cov.val[which(tab.cov.param1[1,1]==parameters)]<-exp(get(tab.cov.param1[get(tab.cov.param1[1,2])==tab.cov.param1[,4],3]))
                      
                  } else 
                  {
                  thetah[which(tab.cov.param1[1,1]==parameters)]<-thetah[which(tab.cov.param1[1,1]==parameters)]+ get(tab.cov.param1[get(tab.cov.param1[1,2])==tab.cov.param1[,4],3])  

                  }
                 }
               } else{  #cas d'un parametre associé à plusieurs covariables
                
                       p_val<-1
                       som<-0
                       for (i.par in 1:length(unique(tab.cov.param1[,2]))) #boucle sur les covaraibles pour le parametres
                       {
                             
                             par.test1<-unique(tab.cov.param1[,2])[i.par]
                             tab.cov.param2<-tab.cov.param1[which(tab.cov.param1[,2]==par.test1),]
                             tab.cov.param2<-matrix(tab.cov.param2,ncol=5)
                             
                             test.bet<-tab.cov.param2[get(tab.cov.param2[1,2])==tab.cov.param2[,4],3]
                           
                             if (length(test.bet)!=0)  #si categories ==0
                            {
                               if (Trand==2)   
                              {
                                 
                                  thetah[which(tab.cov.param2[1,1]==parameters)]<-thetah[which(tab.cov.param2[1,1]==parameters)]*exp(get(tab.cov.param2[get(tab.cov.param2[1,2])==tab.cov.param2[,4],3]))
                                  parameters.cov.val[which(tab.cov.param2[1,1]==parameters)]<-exp(get(tab.cov.param2[get(tab.cov.param2[1,2])==tab.cov.param2[,4],3]))
                                  
                              } else {
                              
                                  thetah[which(tab.cov.param[i,1]==parameters)]<-thetah[which(tab.cov.param[i,1]==parameters)]+ get(tab.cov.param2[get(tab.cov.param2[1,2])==tab.cov.param2[,4],3])  

                              }
                            }
                        }
                   }

                }
                    
            }
        		
  
	for (ip in 1:p) {assign(parameters[ip],thetah[ip])} 
	#for (i in 1:p) {assign(parameters[i],theta[i])}

 
  for (l in 1:nr){
	   s<-matrix(nrow=p,ncol=lprotk[[l]]) 
	   s3<-matrix(nrow=p,ncol=lprotk[[l]])
     
     mod<-numeric() 
      for (i in 1:p){
	         for (j in 1:lprotk[[l]])
		       {   
                #for (ip in 1:p) {assign(parameters[ip],theta_init[ip])} 
                vlf<-1:lformg[[l]]
                forma<-formg[[l]]
                lf<-length(forma)
                Dforma<-lapply(1:lf,function(p,forma,parameters,lf) lapply(1:p,function(p,forma,parameters) D(forma,parameters[p]), forma=forma[lf],parameters=parameters),p=p,forma=forma,parameters=parameters)
               #c<-c(0,cumsum(lformg))
                t<-protk[[l]][j]
                ff1<-formg[[l]][t<=tf[[l]]][1]
			          u<-eval(Dforma[[vlf[t<=tf[[l]]][1]]][[i]])#*dose #
			          #IF COVARIATE  only
			          if(!is.null(categories.cov)) s3[i,j]<-u*parameters.cov.val[i]
                
                if (Trand==2)  u<-beta[i]*u*exp(beta.occh[i])
		            mod[j]<-eval(ff1)#*dose #
		            #IF COVARIATE  only
                if(!is.null(categories.cov)) s[i,j]<-u*parameters.cov.val[i] else s[i,j]<-u
           
		        }
		     
	     }                                       
	   mod1<-c(mod1,mod)
	   s1<-cbind(s1,s)  # Trand==2
	   
	   #IF COVARIATE ONLY
     if(!is.null(categories.cov)) s4<-cbind(s4,s3) # pas de multiplication par valeur du parametre Trand==1
    
  }
  
  #IF COVARIATE ONLY
  if(!is.null(categories.cov)) {
  #covariate parameters
     for (icat in 1: dim(tab.cov.param)[1])
     {
      
      if (!tab.cov.param[icat,4]==get(tab.cov.param[icat,2])) s2<-s1[which(parameters==tab.cov.param[[icat]]),]*0  else   s2<-s1[which(parameters==tab.cov.param[[icat]]),]
       
       
       s1<-rbind(s1,s2)
       s4<-rbind(s4,s2)
     }
      s1rand<-as.matrix(s1[1:length(theta),])
      l1<-list(mod1,s1,s1rand,n.theta,s4)
      thetah<-theta_init
  }else{  
    # (no covariate)                          
    if (dim(s1)[2]==0) l1<-list(mod1,s1,s1,theta_init,s1)
    if (dim(s1)[2]!=0) l1<-list(mod1,s1,s1,theta_init,s1/beta)
  }
	return(l1)

}



#Dériver par rapport aux effets fixes
sensifixe<-function(nr,protk,pfixe,p,kjk,lprotk,doses,categories.cov,names.cov,beta.covariate,theta2,tab.cov.param,n_occ,sequ)
{
  mod<-numeric()
  c<-sum(unlist(lprotk))
  s<-matrix(nrow=pfixe+locc,ncol=0)
  
  for(occ in 1:n_occ){
    
    sensibili<-sensibilite(nr,protk,pfixe,p,kjk,lprotk,doses,categories.cov,names.cov,beta.covariate,theta2[[sequ]][[occ]],tab.cov.param,beta.occ2[[sequ]][[occ]])
    if (Trand ==2) 
    {
      sensibilif<-sensibili[[5]]
      }else 
    { 
      sensibilif<-sensibili[[2]] 
    }
    
    mod_int<-sensibili[[1]]
    mod<-c(mod,mod_int)
    
    bloc0<-matrix(rep(0,locc*c),nrow=locc,ncol=c)
    if (lP[[sequ]][occ]==0) s_int<-rbind(sensibilif,bloc0)

    if (lP[[sequ]][occ]>0){ 
      block<-matrix(nrow=0,ncol=c)
      for (kj in 1:nbcov_occ){
        Seq<-comb.seq.poss[sequ,kj]
        l<-which(covariate_occ.sequence[[kj]][[Seq]][occ]==covariate_occ.category[[kj]])-1 
        block_int<-matrix(0,nrow=length(parameter_occ.associated[[kj]])*(length(covariate_occ.category[[kj]])-1),ncol=c)
        if (length(which(vecP1[[kj]][[Seq]][[occ]]==T))>0) 
        {
          for (j in 1:length(which(vecP1[[kj]][[Seq]][[occ]]==T))){
          i=l+(j-1)*(length(covariate_occ.category[[kj]])-1)
          if (Trand==2) s_int_ij<-as.matrix(sensibili[[2]][1:p,])[which(vecP1[[kj]][[Seq]][[occ]]==T)[j],]
          else  s_int_ij<-as.matrix(sensibili[[5]][1:p,])[which(vecP1[[kj]][[Seq]][[occ]]==T)[j],]
          block_int[i,]<-s_int_ij
        }
      }
      block<-rbind(block,block_int)
       if (Trand ==2 && covariate_occ.model==T) 
    {
    sensibilif[which(vecP1[[kj]][[Seq]][[occ]]==T),]<-sensibilif[which(vecP1[[kj]][[Seq]][[occ]]==T),]*exp(beta.occ2[[sequ]][[occ]][which(vecP1[[kj]][[Seq]][[occ]]==T)])
   
    }
      
    }
    s_int<-rbind(sensibilif,block)
   }                        
   s<-cbind(s,s_int)
  }
  l<-list(mod,s)
  return(l)
}


#Dériver par rapport aux effets aléatoires
sensialea<-function(nr,protk,pfixe,p,kjk,lprotk,doses,categories.cov,names.cov,beta.covariate,theta2,tab.cov.param,n_occ,sequ)
{
  c<-sum(unlist(lprotk))
  mod<-numeric()
  s1<-matrix(nrow=p,ncol=0)
  s2<-matrix(0,nrow=p*pn,ncol=c*pn)
  s<-matrix()
   
  if (n_occ==1) s<-sensibilite(nr,protk,pfixe,p,kjk,lprotk,doses,categories.cov,names.cov,beta.covariate,theta2[[sequ]][[1]],tab.cov.param,beta.occ2[[sequ]][[1]])[[3]]
  else{ 
    for(occ in 1:pn)
    { sensibili<-sensibilite(nr,protk,pfixe,p,kjk,lprotk,doses,categories.cov,names.cov,beta.covariate,theta2[[sequ]][[occ]],tab.cov.param,beta.occ2[[sequ]][[occ]])
      mod_int<-sensibili[[1]]
      s_int<-sensibili[[3]]
      mod<-c(mod,mod_int)
      s1<-cbind(s1,s_int)
      s2[(p*(occ-1)+1) : (p*occ),(c*(occ-1)+1) :(c*occ)]=s_int
    }
    s<-rbind(s1,s2)
  }
  l<-list(mod,s)
  return(l)
}


# Fisher matrix for 1 ind in 1 elementary design by block: block A (resA), block B(resB), block C (resC)
# avec IOV
#fonction pour écrire une matrice block diagonal à partir d'une liste de matrices
blockdiag <-function(x){
     if(!is.list(x)) stop("x not a list")
     n <- length(x)
     if(n==0) return(NULL)
     x <- lapply(x, function(y) if(length(y)) as.matrix(y) else
stop("Zero-length component in x"))
     d <- array(unlist(lapply(x, dim)), c(2, n))
     rr <- d[1,]
     cc <- d[2,]
     rsum <- sum(rr)
     csum <- sum(cc)
     out <- array(0, c(rsum, csum))
     ind <- array(0, c(4, n))
     rcum <- cumsum(rr)
     ccum <- cumsum(cc)
     ind[1,-1] <- rcum[-n]
     ind[2,] <- rcum
     ind[3,-1] <- ccum[-n]
     ind[4,] <- ccum
     imat <- array(1:(rsum * csum), c(rsum, csum))
     iuse <- apply(ind, 2, function(y, imat) imat[(y[1]+1):y[2],
(y[3]+1):y[4]], imat=imat)
     iuse <- as.vector(unlist(iuse))
     out[iuse] <- unlist(x)
     return(out)
} 

#-------------------------------------------------------------------------------------------------------
bloc<-function(protk,sensif,sensia,var,invvar,pfixe,p,dvar,mod,lprotk,bout4,categories.cov,sequ)   #proti,protiPD,sensi,var,invvar,dvar,mod
{

#Block A
#-------
  resA<-matrix(nrow=pfixe+locc,ncol=pfixe+locc)
  resA<-(sensif %*% invvar) %*% t(sensif)  #matrice parametres fixe#

  #resA<-(cbind(sensibili[[5]],sensibili[[5]]) %*% invvar) %*% t(cbind(sensibili[[5]],sensibili[[5]]))
#Block B
#-------
  n_occB<-n_occ
  if (length(which(gamma==0))==length(gamma)) n_occB<-1
  if (n_occB>1) l<-2*p else l<-p
	resB<-matrix(nrow=l+nr*2,ncol=l+nr*2)     #matrice aléatoire (block B) #
	resB_int<-matrix(nrow=l+nr*2,ncol=l+nr*2)     #calcul intermédiaire
	#resC<-matrix(rep(0,(p+nr*2)*pfixe),nrow=p+nr*2,ncol=pfixe)     #matrice incomplete option 1 covariance==0#

	#if (Trand==2)  sensia<-sensi.theta else sensia<-sensi 
  
	#resA1<-(sensia %*% invvar) %*% t(sensia)  #matrice parametres fixe#
	
	l1<-list()
	li<-list()
	l2<-list()
	
	#premiere réponse
	c<-list()
	long<-c()
	a<-NULL
	for (occ in 1:n_occB){
	c[[occ]]<-cumsum(unlist(lprotk))+length(unlist(protk))*(occ-1)
	long[occ]<-c[[occ]][nr]-c[[occ]][1]
	a0<-c(bout4[[occ]][[1]],rep(0,long[occ]))
	a<-c(a,a0)
	}
	a1<-2*diag(a,length(unlist(protk))*n_occB)
	ai<-invvar%*%a1%*%invvar
	resB_int[l+1,l+1]<-sum(diag(ai%*%a1))

	a2<-a1*mod
	resB_int[l+2,l+2]<-sum(diag(a2%*%invvar%*%a2%*%invvar))
	resB_int[l+1,l+2]<-sum(diag(ai%*%a2))
	resB_int[l+2,l+1]<-resB_int[l+1,l+2]
	l1[[1]]<-a1
	li[[1]]<-ai
	l2[[1]]<-a2
	
	
	if (nr!=1){
	for (i in 1:(nr-1)){
	vec<-NULL
	 for (occ in 1:n_occB){
	   vec0<-rep(0,length(unlist(protk)))
	   if(length(protk[[i+1]])!=0){
	   
	     vec0<-c(rep(0,length(protk[[i]])),bout4[[occ]][[i+1]])
	   } else {vec0<-vec0}
	 
   vec<-c(vec,vec0)
   }
   
	a1<-2*diag(vec,length(unlist(protk))*n_occB)

	ai<-invvar%*%a1%*%invvar
	resB_int[l+i*2+1,l+i*2+1]<-sum(diag(ai%*%a1))
	#sigma slope #
	a2<-a1*mod
  resB_int[l+i*2+2,l+i*2+2]<-sum(diag(a2%*%invvar%*%a2%*%invvar))
	resB_int[l+(i*2)+1,l+i*2+2]<-sum(diag(ai%*%a2))
	resB_int[l+i*2+2,l+(i*2)+1]<-resB_int[l+(i*2)+1,l+i*2+2]
	l1[[i+1]]<-a1
	li[[i+1]]<-ai
	l2[[i+1]]<-a2
  }
 
 
  #sigma inter PD et sigma inter PK#
  for (i in 1: (nr)){
    for (j in  2 : (nr)){

	     resB_int[(l+j*2-1),l+i*2-1]<-sum(diag(li[[j]]%*%l1[[i]])) 
       resB_int[l+i*2-1,(l+j*2-1)]<-resB_int[(l+j*2-1),l+i*2-1] #resB[p+1,p+3]<-resB[p+3,p+1]
	     resB_int[l+j*2,l+i*2-1]<-sum(diag(li[[i]]%*%l2[[j]])) 	#resB[p+1,p+4]<-resB[p+4,p+1]
       resB_int[l+i*2-1,l+j*2]<-resB_int[l+j*2,l+i*2-1]
	    
       resB_int[l+j*2-1,l+i*2]<-sum(diag(li[[j]]%*%l2[[i]]))  #resB[p+2,p+3]<-resB[p+3,p+2]
	     resB_int[l+i*2,l+j*2-1]<-resB_int[l+j*2-1,l+i*2]
       resB_int[l+j*2,l+i*2]<-sum(diag(invvar%*%l2[[j]]%*%invvar%*%l2[[i]]))
       resB_int[l+i*2,l+j*2]<-resB_int[l+j*2,l+i*2] #resB[p+2,p+4]<-resB[p+4,p+2]
	}
	}
	
	for (i in 1:p)
	{
    resB1<-(sensia[i,])%*%t(sensia[i,])
    resB2<-invvar%*%resB1%*%invvar
    resB_int[i,1:p]<-sapply(1:p,function(resB2,sensia,p) sum(diag(resB2%*%(sensia[p,])%*%t(sensia[p,]))),resB2=resB2,sensia=sensia)
    resB_int[i,l+1]<-sum(diag(l1[[1]]%*%resB2))
		resB_int[i,l+2]<-sum(diag(l2[[1]]%*%resB2))
    for (j in 2:nr) {
		resB_int[i,l+j*2-1]<-sum(diag(l1[[j]]%*%resB2))
		resB_int[i,l+j*2]<-sum(diag(l2[[j]]%*%resB2)) 	
    }
		resB_int[(l+1):(l+(2*nr)),i]<-resB_int[i,(l+1):(l+(2*nr))]
		
    if (n_occB>1 && length(which(gamma==0))!=length(gamma))
    {
    resB3<-blockdiag(sapply(1:n_occB,function(mat,n_occB) list(((sensia[i,])%*%t(sensia[i,]))[1:length(unlist(protk)),1:length(unlist(protk))]),mat=((sensia[i,])%*%t(sensia[i,]))[1:length(unlist(protk)),1:length(unlist(protk))]))
    resB4<-invvar%*%resB3%*%invvar
    resB_int[p+i,1:p]<-sapply(1:p,function(resB4,sensia,p) sum(diag(resB4%*%(sensia[p,])%*%t(sensia[p,]))),resB4=resB4,sensia=sensia)
    resB_int[1:p,p+i]<-resB_int[p+i,1:p]
    resB_int[p+i,(p+1):l]<-sapply(1:p,function(resB4,sensia,p) sum(diag(resB4%*%blockdiag(sapply(1:n_occB,function(mat,n_occB) list(((sensia[p,])%*%t(sensia[p,]))[1:length(unlist(protk)),1:length(unlist(protk))]),mat=((sensia[p,])%*%t(sensia[p,]))[1:length(unlist(protk)),1:length(unlist(protk))])))),resB4=resB4,sensia=sensia)
		resB_int[p+i,l+1]<-sum(diag(l1[[1]]%*%resB4))
		resB_int[p+i,l+2]<-sum(diag(l2[[1]]%*%resB4))
    for (j in 2:nr) {
		resB_int[p+i,l+j*2-1]<-sum(diag(l1[[j]]%*%resB4))
		resB_int[p+i,l+j*2]<-sum(diag(l2[[j]]%*%resB4)) 	
    }
		resB_int[(l+1):(l+(2*nr)),p+i]<-resB_int[p+i,(l+1):(l+(2*nr))]
		}
	}
# Déduire B à partir de B_intermédiare
 resB<-resB_int 	
 }
 

  
  else {
	 for (i in 1:p)
	{
    resB1<-(sensia[i,])%*%t(sensia[i,])
    resB2<-invvar%*%resB1%*%invvar
    resB_int[i,1:p]<-sapply(1:p,function(resB2,sensia,p) sum(diag(resB2%*%(sensia[p,])%*%t(sensia[p,]))),resB2=resB2,sensia=sensia)
    resB_int[i,l+1]<-sum(diag(l1[[1]]%*%resB2))
		resB_int[i,l+2]<-sum(diag(l2[[1]]%*%resB2))
		resB_int[(l+1):(l+2),i]<-resB_int[i,(l+1):(l+2)]
		
    if (n_occB>1 && length(which(gamma==0))!=length(gamma))
    {
    resB3<-blockdiag(sapply(1:n_occB,function(mat,n_occB) list(((sensia[i,])%*%t(sensia[i,]))[1:length(unlist(protk)),1:length(unlist(protk))]),mat=((sensia[i,])%*%t(sensia[i,]))[1:length(unlist(protk)),1:length(unlist(protk))]))
    resB4<-invvar%*%resB3%*%invvar
    resB_int[p+i,1:p]<-sapply(1:p,function(resB4,sensia,p) sum(diag(resB4%*%(sensia[p,])%*%t(sensia[p,]))),resB4=resB4,sensia=sensia)
    resB_int[1:p,p+i]<-resB_int[p+i,1:p]
    resB_int[p+i,(p+1):l]<-sapply(1:p,function(resB4,sensia,p) sum(diag(resB4%*%blockdiag(sapply(1:n_occB,function(mat,n_occB) list(((sensia[p,])%*%t(sensia[p,]))[1:length(unlist(protk)),1:length(unlist(protk))]),mat=((sensia[p,])%*%t(sensia[p,]))[1:length(unlist(protk)),1:length(unlist(protk))])))),resB4=resB4,sensia=sensia)
    resB_int[p+i,l+1]<-sum(diag(l1[[1]]%*%resB4))
		resB_int[p+i,l+2]<-sum(diag(l2[[1]]%*%resB4))
		resB_int[(l+1):(l+2),p+i]<-resB_int[p+i,(l+1):(l+2)]
		}
	}
  # Déduire B à partir de B_intermédiare
 resB<-resB_int  	   
}

  
  # B final
  resB<-0.5*resB                               
  # Bloc C 
  resC<-matrix(rep(0,(l+nr*2)*(pfixe+locc)),nrow=l+nr*2,ncol=pfixe+locc)     #matrice incomplete option 1 covariance==0#
	# où pfixe est 1 dim de A	
    
  h<-cbind(rbind(resA,resC),rbind(t(resC),resB))
	return(h) 
}


#Population Fisher information matrix avec iov
#------------------------------------
fisher<-function(prots,subjects,doses,Ntot,ind=0,nq,categories,names.cov,proportions.cov,beta.covariate,tab_mod_prop,pfixe,p,theta2,tab.cov.param,n_occ,lposs)
{	 
  n_occB<-n_occ
  if (length(which(gamma==0))==length(gamma)) n_occB<-1
	if (n_occB>1) pp<-pfixe+length(beta.occ1)+p*2+2*nr  else pp<-pfixe+length(beta.occ1)+p+2*nr
	somme<-matrix(c(rep(0,pp*pp)),ncol=pp) 
	   #test if covariate or not
	   #categories <-covariate.cat.poss
#names.cov<-info.covariate[[1]]
  no.cov<-1
  
    if (is.null(categories)) 
        {
        categories<-matrix(0, nrow=1,ncol=1)
		    no.cov<-NULL
		    }
	#For each elementary design computation of the elementary fisher information matrix
for (sequ in 1:lposs){
    #somme1<-matrix(c(rep(0,pp*pp)),ncol=pp) 
	 for  (i in 1:length(prots[[1]]))	#meme nombre de temps pour chaque reponse par groupe elementaire
		{  
		    
		    #categories:table with for each line a combination of the modalities for each covariates
		    #loop on all possible combination of modalities of all covariates 

		    for (ll in 1:dim(categories)[1])
        {
           
           if (is.null(no.cov)) categories.cov<-NULL else categories.cov<-c(categories[ll,])
           if (is.null(no.cov)) categories.prop<-1 else categories.prop<-prod(tab_mod_prop[[2]][ll,])  #si no covariate categories.prop=1

	     	   protk<-pts(nr,i,prots)
           lprotk<-lapply(1:nr,function(nr,protk) length(protk[[nr]]),protk=protk )
           sensi<-sensifixe(nr,protk,pfixe,p,i,lprotk,doses,categories.cov,names.cov,beta.covariate,theta2,tab.cov.param,n_occB,sequ) 
           mod<-sensi[[1]] 
           sensif<-sensi[[2]]   #matrice derivée premiere effets fixes (avec covariables si demandé)
           sensia<-sensialea(nr,protk,pfixe,p,i,lprotk,doses,categories.cov,names.cov,beta.covariate,theta2,tab.cov.param,n_occB,sequ)[[2]]
           if (n_occB>1) bout1<-OMEGA%*%sensia else bout1<-omega%*%sensia
           
			     bout4<-list()
           t1<-0
	      	 t<-0
	      	 for (occ in 1:n_occB){
      	   bout4[[occ]]<-list()
			     for (n in 1:(nr))
           {
		    	     t<-1+t1
			         t1<-t-1+lprotk[[n]]
			
			         if(lprotk[[n]]==0)
                {
			            bout4[[occ]][[n]]<-NULL 
                #if(n==nr){bout4[[occ]][[n+1]]=Inf}
                } 
                else 
                {
		        	     bout3<-sigmainter[[n]]+sigmaslope[[n]]*mod[t:t1]
		        	     bout4[[occ]][[n]]<-bout3
		            }
           }
           }
           
		    	 t<-NULL
			     t1<-NULL
			     var<-t(bout1) %*% sensia + diag(unlist(bout4)^2,length(unlist(protk))*n_occB) 
		        
			     if (length(which(dim(var)==0))!=2) {invvar<-solve(as.matrix(var))
           fish<-bloc(protk,sensif,sensia,var,invvar,pfixe,p,dvar,mod,lprotk,bout4,categories.cov,sequ)}
           else {fish<-matrix(c(rep(0,pp*pp)),ncol=pp)} 
           
           somme<-somme+fish*categories.prop*comb.seq.prop[sequ]*subjects[i]     
        }    
			 }
	   }
 
		sigmaslopen<-lapply(1:nr,function(nr,x) x[nr],x=paste("sig.slope",LETTERS[1:nr],sep=""))
    sigmaintern<-lapply(1:nr,function(nr,x) x[nr],x=paste("sig.inter",LETTERS[1:nr],sep=""))
    xn<-rbind(sigmaintern,sigmaslopen)
    xv<-rbind(sigmainter,sigmaslope)
    SIG<<-unlist(c(xv))
    SIGN<<-unlist(c(xn))
    vec2<<-SIG!=0
	  SIG<<-SIG[vec2]	
	  SIGN<-SIGN[vec2]
	  lSIG<<-length(SIG)		

    OMEGABIS<-diag(c(diag(omega),diag(gamma)))
    OMEGABIS<-as.matrix(OMEGABIS)
    if (length(c(diag(omega),diag(gamma)))==1) OMEGABIS<-as.matrix(c(diag(omega),diag(gamma)))
		if (n_occB>1) vec1<<-diag(OMEGABIS)!=0 
		if (n_occB==1) vec1<<-diag(omega)!=0 
    OMEGA1<<-as.matrix(as.matrix(OMEGABIS[,vec1])[vec1,]) 
   	omega1<<-as.matrix(as.matrix(omega[,diag(omega)!=0])[diag(omega)!=0,])
   	gamma1<<-as.matrix(as.matrix(gamma[,diag(gamma)!=0])[diag(gamma)!=0,])
    if (n_occB==1) OMEGA1<-omega1
    
		vec<-c(rep(T,pfixe+locc),vec1,vec2) #remove the row i and the col i of somme if omega[i,i]=0  
		somme<-somme[vec,vec]

    pp<-dim(somme)[1] 	
    det<-det(somme)	 
		inv<-try(solve(as.matrix(somme))) 
	   if(!is.null(attributes(inv)$class))
	    {
		se<-rep(NA,pp)
		cv<-se	
		crit<-(det)^(1/pp) 
		} else
     {
		se<-sqrt(diag(inv))   
    if (n_occ>1) {  
		cv1<-se[1:pfixe]/c(theta,unlist(beta.covariate))*100 
    cv1bis<-se[(pfixe+1):(pfixe+locc)]/beta.occ1*100 
		cv2<-se[(pfixe+locc+1):(pp-lSIG)]/diag(as.matrix(OMEGA1))*100 
		cv3<-se[(pp-(lSIG-1)):pp]/SIG*100 
		cv<-abs(c(cv1,cv1bis,cv2,cv3))	
		}
		else{
		cv1<-se[1:pfixe]/c(theta,unlist(beta.covariate))*100 
		cv2<-se[(pfixe+1):(pp-lSIG)]/diag(as.matrix(OMEGA1))*100
		#if (n_occ==1 && nr==1 && dim(!matrix(OMEGA1))==c(1,1))   cv2<-se[(pfixe+1):(pp-lSIG)]/c(OMEGA1)*100 
		cv3<-se[(pp-(lSIG-1)):pp]/SIG*100 
		cv<-abs(c(cv1,cv2,cv3))	
    }
		crit<-(det)^(1/pp) 
	   }
		return(list(form,somme,subjects,doses,prots,se,cv,det,crit,pfixe,pp,xn,p)) 
}


#Calculer l'intervalle de confiance
###################################

CI.computation<-function(alpha,beta.cov,beta.cov.se,nom.par,Trand)
{
norm.val<-abs(qnorm(alpha, mean=0, sd=1))

Beta<-beta.cov
lcov<-length(beta.cov)
#confidence interval for beta
CI_inf<-beta.cov-norm.val*beta.cov.se
CI_sup<-beta.cov+norm.val*beta.cov.se 
t<-unlist(lapply(1:lcov,function(lcov) paste("[",round(CI_inf[lcov],3),";",round(CI_sup[lcov],3),"]",sep=""))) 
lim<-paste(100*(1-2*alpha),"% CI") 
if (Trand==2)
{
#Confidence interval for exp(beta)   
Beta_exp<-exp(beta.cov)
vide<-rep(" ",lcov) 
t1<-unlist(lapply(1:lcov,function(lcov) paste("[",round(exp(CI_inf[lcov]),3),";",round(exp(CI_sup[lcov]),3),"]",sep=""))) 
ci<-data.frame(Beta,t,vide,Beta_exp,t1,row.names=paste(nom.par))
names(ci)<-c("Beta",lim," ","exp(Beta)",lim)
}
else
{ci<-data.frame(Beta,t,row.names=paste(nom.par))
names(ci)<-c("Beta",lim)
}

return(ci)

}


#Calculer la puissance et le nombre de sujet pour le test de comparaison
power_nni.computation<-function(alpha,given.power,init.num.subjects,beta.cov,beta.cov.se,nom.par,Trand)
{
norm.val<-abs(qnorm(alpha, mean=0, sd=1))

#Puissance avec cet effet la avec alpha =5%
power.prediction<-1-pnorm(norm.val-(beta.cov/beta.cov.se)) + pnorm(-norm.val-(beta.cov/beta.cov.se))

Beta<-beta.cov
Predicted_SE<-beta.cov.se
lcov<-length(beta.cov)  
Expected_power<-round(power.prediction,6)

#Calcul du nombre de sujet necessaires pour avoir une puissance donnée%
#given.power
#calcul de la SE si je me met sous H1
NSE<-beta.cov/(norm.val-qnorm(1-given.power))

#calcul maintenant du nombre de sujet necessaires pour mettre en evidence cette effets
#avec une puissance donnée
NNI<-sum(init.num.subjects)*(beta.cov.se/NSE)^2

Number_subjects_needed<-NNI

#Affichage
vide<-rep(" ",lcov)
if (compute.power==T){
result<-data.frame(Expected_power,row.names=paste(nom.par))
names(result)<-c("Expected_power")
if (compute.nni==T){
result<-cbind(result,Number_subjects_needed,vide)
names(result)<-c("Expected_power","Number_subjects_needed",paste("(for a given power=",given.power,")",sep=""))
}
}
if (compute.power==F && compute.nni==T){
result<-data.frame(Number_subjects_needed,vide,row.names=paste(nom.par))
names(result)<-c("Number_subjects_needed",paste("(for a given power=",given.power,")",sep=""))
}
return(result)
}



#Calculer la puissance et le NSN pour le test d'équivalence
power_nni_eq.computation<-function(alpha,given.power,init.num.subjects,beta.cov,beta.cov.se,nom.par,Trand)
{
norm.val<-abs(qnorm(alpha, mean=0, sd=1))
Beta<-beta.cov
Predicted_SE<-beta.cov.se
lcov<-length(beta.cov)

#puissance
power.prediction<-c()
for (i in 1:length(beta.cov)){
if(beta.cov[i]==0){power.prediction[i]<-pnorm(-norm.val+(interval_eq[2]/beta.cov.se[i]))}
if(beta.cov[i] <0){power.prediction[i]<-1-pnorm(norm.val-(beta.cov[i]-interval_eq[1])/beta.cov.se[i])}
if(beta.cov[i] >0){power.prediction[i]<-pnorm(-norm.val-(beta.cov[i]+interval_eq[1])/beta.cov.se[i])}
#if (beta.cov[i] < interval_eq[1] || beta.cov[i] > interval_eq[2]) {power.prediction[i]<-c(".")
#test_eq<-F}
}




#Calcul du nombre de sujet necessaires pour avoir une puissance donnée%
#given.power
#calcul de la SE si je me met sous H1
NSE<-c()
for (i in 1:length(beta.cov)){
if (beta.cov[i]==0){NSE[i]<-interval_eq[2]/(norm.val+qnorm(given.power))}
if (beta.cov[i]<0){NSE[i]<-(-beta.cov[i]+interval_eq[1])/(-norm.val+qnorm(1-given.power))}
if (beta.cov[i]>0){NSE[i]<-(-beta.cov[i]-interval_eq[1])/(norm.val+qnorm(given.power))}
}
#calcul maintenant du nombre de sujet necessaires pour mettre en evidence cette effets
#avec une puissance donnée
NNI<-sum(init.num.subjects)*(beta.cov.se/NSE)^2



#si pas sous H1

for (i in 1:length(beta.cov)){
if (beta.cov[i] < interval_eq[1] || beta.cov[i] > interval_eq[2]) {
power.prediction[i]<-c(NA)
NNI[i]<-c(NA)
}#else
#{#if (length(which(beta.cov[i] < interval_eq[1])) + length(which(beta.cov[i]> interval_eq[2])) < length(beta.cov))
#power.prediction[i]<-round(power.prediction[i], digits = 6)
#NNI[i]<-round(NNI[i], digits = 6)} 
}

Expected_power<-power.prediction
Number_subjects_needed<-NNI
#Affichage
vide<-rep(" ",lcov)

if (compute.power_eq==T){
result<-data.frame(Expected_power,row.names=paste(nom.par))
names(result)<-c("Expected_power")
if (compute.nni_eq==T){
result<-cbind(result,Number_subjects_needed,vide)
names(result)<-c("Expected_power","Number_subjects_needed",paste("(for a given power=",given.power,")",sep=""))
}
}
if (compute.power_eq==F && compute.nni_eq==T){
result<-data.frame(Number_subjects_needed,vide,row.names=paste(nom.par))
names(result)<-c("Number_subjects_needed",paste("(for a given power=",given.power,")",sep=""))
}
return(result)
}


#GRAPHES

graph<-function(formg,prots,graphinf,graphsup,names.datax,names.datay,y.range)
{

		for (i in 1:p) {assign(parameters[i],theta[i]) }

		if (dose.identical==T)
		{
		    dose<-doses[1] 
		    ff1<-list()
		    ff2<-list()
		    
        par(mfrow=c(1,nr))
        for (l in 1 : nr){
            
		        ff1[[l]]<-vector(mode="numeric", length(tt1[[l]]))
		        
				    for (j in 1:length(tt1[[l]]))
	          {
              #ff1[[l]]<-vector(mode="numeric",length=length(tt[[l]]))
				      t<-tt1[[l]][j]
				         ff1[[l]][j]<-eval(formg[[l]][t<=tf[[l]]][1])
    	      
  		      }
  		      
  		      
  		      if(log.logical!=F) plot(tt1[[l]],ff1[[l]],type='l',xlab=paste(names.datax[l]),ylab=paste(names.datay[l]),log=log.logical,ylim=c(y.range[[l]],range(ff1[[l]]))[1:2]) 
            else   plot(tt1[[l]],ff1[[l]],type='l',xlab=paste(names.datax[l]),ylab=paste(names.datay[l]),ylim=c(y.range[[l]],range(ff1[[l]]))[1:2])
        
            #points for the design on the graph
				          
				          for ( i in 1:length(prots[[l]]))
				          {
				          ff2[[l]]<-vector(mode="numeric", length(prots[[l]][[i]]))
				            for ( j in 1:length(prots[[l]][[i]])) {
				    	          t<-prots[[l]][[i]][j]
					              ff2[[l]][j]<-eval(formg[[l]][t<=tf[[l]]][1])#*dose[1] #
		                }
				                points(prots[[l]][[i]],ff2[[l]],pch=paste(i),cex=2,col=paste(i))
	                }
             
		      }
      }
		else
		{
	    ff1<-list()
	    par(mfrow=c(nr,length(prots[[1]])))
			for (l in 1 : nr){

			 ff1[[l]]<-vector(mode="numeric", length(tt1[[l]]))
			     for (kjk in 1:length(prots[[l]]))
			     {
			#windows()
			#par(mfrow=c(1,2))
			       dose<-doses[kjk]
				      for (j in 1:length(tt1[[l]]))
				      {
				      	 t<-tt1[[l]][j]
					       ff1[[l]][j]<-eval(formg[[l]][t<=tf[[l]]][1])#*dose[k]  #
                
				      }
	         
				      if (log.logical!=F) 
				          plot(tt1[[l]],ff1[[l]],type='l',xlab=paste(names.datax[l]),ylab=paste(names.datay[l]),log=log.logical,ylim=c(y.range[[l]],range(ff1[[l]]))[1:2]) else plot(tt1[[l]],ff1[[l]],type='l',xlab=paste(names.datax[l]),ylab=paste(names.datay[l]),ylim=c(y.range[[l]],range(ff1[[l]]))[1:2])
				
        
              ff2<-list()
              ff2[[l]]<-vector(mode="numeric", length(prots[[l]][[kjk]]))
				      for ( i in 1:length(prots[[l]][[kjk]]))
				      {
	               dose<-doses[kjk]
					       t<-prots[[l]][[kjk]][i]
					       ff2 [[l]][i]<-eval(formg[[l]][t<=tf[[l]]][1])#*dose[k] #				
		          }
           	    points(prots[[l]][[kjk]],ff2[[l]],pch=paste(kjk),cex=2)
			         	title(paste("Group",kjk,"- dose = ", doses[kjk],sep=" "))
			}
			}
      	
		}		
}


out<-function()
{
  d<-Sys.time()  
  if (pn==1 || locc==0) no<-NULL
      if (n_occ>1 && locc>0){
      no<-NULL
      no1<-NULL
      lo<-NULL
  if (covariate_occ.model==T){
  
      #TEST COVARIATE_OCC MODEL
          if ((length(covariate_occ.name)!= length(covariate_occ.category))==T) 
          stop(" The list of the name of the covariate(s) depending also on the occasion has not the same length of the list of the category")
          
          if(length(which((sapply(covariate_occ.sequence,length) == sapply(sequence.proportions,length))==F))!=0 )
          stop(" The number of covariate categorie sequences is not equal to the number of proportions")
          
          if (length(covariate_occ.name)!= length(parameter_occ.associated))
          stop(" One or more covariate(s) has not at least one associated parameter ")
          
          if(sum(unlist(sapply(sequence.proportions,sum)))!=length(covariate_occ.name))
          stop(" The sum of the proportions of the categories is not equal to one")
      

      count_occ<-0
      for (kj in 1:length(covariate_occ.name))
      {
          for (j in 1:length(parameter_occ.associated[[kj]]))
          {
            
              lo<-length(covariate_occ.category[[kj]])-1
              no1<-lapply(1:lo,function(lo) paste("beta_",parameter_occ.associated[[kj]][j],"_",covariate_occ.name[[kj]],"_",covariate_occ.category[[kj]][lo+1],sep=""))
              no<-c(no,unlist(no1)) #liste des parametres
            }
             count_occ<-count_occ+lo
          }

      }
   }




  if(covariate.model==T)
  { #TEST COVARIATE MODEL
          if ((length(covariate.name)!= length(covariate.category))==T || (length(covariate.category)!=length(covariate.proportions))==T) 
          stop(" The list of the name of the covariate(s) has not the same length of the list of the category or the proportions")
          
          if(length(which((sapply(covariate.category,length) == sapply(covariate.proportions,length))==F))!=0 )
          stop(" The number of categories is not equal to the number of proportions")
          
          if (length(covariate.name)!= length(parameter.associated))
          stop(" One or more covariate(s) has not at least one associated parameter ")
          
          if(sum(unlist(sapply(covariate.proportions,sum)))!=length(covariate.name))
          stop(" The sum of the proportions of the categories is not equal to one")
          
        
          covariate.cat.reco<-list()
          info.covariate<-list(covariate.name,covariate.category,covariate.proportions,parameter.associated,beta.covariate)
          info.power.nni<-list(alpha,compute.power,compute.nni,given.power)
          
          #recodage covariable si en character en vecteur numerique#METTRE UN lAPPLY ???
          for (i in 1: length(info.covariate[[1]]))
          {
          if(is.character(info.covariate[[2]][[i]]))
             covariate.cat.reco[[i]]<-c(1:length(info.covariate[[2]][[i]]))
          else 
             covariate.cat.reco[[i]]<-info.covariate[[2]][[i]]
          }
          
          #solution possibles pour les covariables
          covariate.cat.poss<-unique(t(combn(unlist(covariate.cat.reco),length(info.covariate[[1]]))))
          #verifier que les modalités presentes dans chque colonne sont bien des modalités possibles de la covariables
          #test si les modalités possibles pour cahque covaraibles ne sont pas abbérantes
          ll.cov<-length(info.covariate[[1]])
          for (i in 1: ll.cov)
          {
            #nline<-which(covariate.cat.poss[,i] > max(info.covariate[[2]][[i]]))
            nline<-which(covariate.cat.poss[,i] > max(covariate.cat.reco[[i]]) | covariate.cat.poss[,i] < min(covariate.cat.reco[[i]]))
            ref<-integer(length = 0)
        
            if (!identical(nline,ref)) covariate.cat.poss<-covariate.cat.poss[-nline,]
          }
          
          #creation of a table with the parameters aand the covaraite associated 
          npar<-NULL
          npar1<-NULL
          lmodal<-NULL
          count<-0
          ilnp<-0
          s.cov<-matrix(nrow=length(unlist(info.covariate[[5]])),ncol=2)
          for (ijk in 1:length(info.covariate[[1]]))
          {
            for (ijk1 in 1:length(info.covariate[[4]][[ijk]]))
            {
              
              #s.cov[count+ijk1,2]<-info.covariate[[1]][[ijk]]
              #s.cov[count+ijk1,1]<-info.covariate[[4]][[ijk]][ijk1]
              s.cov[count+1,2]<-info.covariate[[1]][[ijk]]
              s.cov[count+1,1]<-info.covariate[[4]][[ijk]][ijk1]
              lnp<-length(info.covariate[[2]][[ijk]])-1
              npar1<-lapply(1:lnp,function(lnp) paste("beta_",info.covariate[[4]][[ijk]][ijk1],"_",info.covariate[[1]][[ijk]],"_",info.covariate[[2]][[ijk]][lnp+1],sep=""))
              npar<-c(npar,unlist(npar1)) #liste des parametres
              
               #liste des modalités
              lmodal1<-lapply(1:lnp,function(lnp) c(covariate.cat.reco[[ijk]][lnp]+1))
              lmodal<-c(lmodal,unlist(lmodal1))
              
              
              if(lnp>1)
              {
                for (i in 2:(lnp))
                { 
                  ilnp<-ilnp+1
                  s.cov[count+1+ilnp,2]<-info.covariate[[1]][[ijk]]
                  s.cov[count+1+ilnp,1]<-info.covariate[[4]][[ijk]][ijk1] 
                  ilnp<-0
                }
              }
             count<-count+lnp
            }
            #count<-count+length(info.covariate[[4]][[ijk]])+ilnp 
            #ilnp<-0
          }
            tab.cov.param<-cbind(s.cov,npar,lmodal,unlist(info.covariate[[5]]))
            
            ###tranform in matrix
            #tab.cov.param1<-tab.cov.param[order(tab.cov.param[,1]),]
            #tab.cov.param<-matrix(tab.cov.param1,nrow=dim(tab.cov.param)[1],ncol=dim(tab.cov.param)[2])
            
            tab.cov.param<-matrix(tab.cov.param,nrow=dim(tab.cov.param)[1],ncol=dim(tab.cov.param)[2])
            
            #number of fixed effects with the covariate
            pfixe<-length(c(parameters,npar))
        
            #summary info covariate model
            #association des prportions aux différentes combianaisons
            tab_mod_prop<-list()
            tab_mod_prop<-association_modality_covariate(info.covariate[[1]],covariate.cat.reco,info.covariate[[3]],covariate.cat.poss)
            f<-fisher(prots,subjects,doses,Ntot,ind=0,nq,covariate.cat.poss,info.covariate[[1]],info.covariate[[3]],info.covariate[[5]],tab_mod_prop,pfixe,p,theta2,tab.cov.param,n_occ,lposs) 
      #fisher2<-function(prots,subjects,doses,Ntot,ind=0,nq,categories,names.cov,proportions.cov,beta.covariate,tab_mod_prop,pfixe,p,theta,theta2,tab.cov.param,pn)

      } else{
           covariate.cat.poss<-NULL
           info.covariate<-list(NULL,NULL,NULL,NULL,NULL)
           tab_mod_prop<-NULL
           pfixe<-length(parameters)
           tab.cov.param<-NULL
          #fisher 2 prend en compte iov 
    f<-fisher(prots,subjects,doses,Ntot,ind=0,nq,covariate.cat.poss,info.covariate[[1]],info.covariate[[3]],info.covariate[[5]],tab_mod_prop,pfixe,p,theta2,tab.cov.param,n_occ,lposs) 
     
      }

  subjects<-f[[3]]
	pfixe<-f[[10]]
	pp<-f[[11]]
	se<-f[[6]] 
	cv<-f[[7]]
	mfisher<-f[[2]]
	determinant<-f[[8]] 
	crit<-f[[9]] 
	StdError<-se[1:(pfixe+locc)] 
	RSE<-cv[1:(pfixe+locc)] 
	
	if (covariate.model==T) 
  {
  Beta<-c(theta,unlist(beta.covariate),beta.occ1) 
  if (n_occ==0) a<-data.frame(Beta,StdError,RSE,row.names=c(parameters,npar)) 
  if (n_occ>0) a<-data.frame(Beta,StdError,RSE,row.names=c(parameters,npar,no))  
  }else{
  Beta<-c(theta,beta.occ1) 
  if (n_occ==0) a<-data.frame(Beta,StdError,RSE,row.names=c(parameters)) 
  if (n_occ>0) a<-data.frame(Beta,StdError,RSE,row.names=c(parameters,no)) 
	}
  .<-c(rep("%",pfixe+locc)) 
	a<-cbind(a,.)
  names(a)[dim(a)[2]]<-" " 
	
	
 
	b_condition<-F
	if (length(which(omega==0))!=length(omega)) {
	Omega<-diag(as.matrix(omega1)) 
	StdError<-se[(pfixe+locc+1):(pfixe+locc+length(Omega))] 
	RSE<-cv[(pfixe+locc+1):(pfixe+locc+length(Omega))]
	b<-data.frame(Omega,StdError,RSE,row.names=parameters[diag(omega)!=0]) 
	.<-c(rep("%",(length(Omega)))) 
	b<-cbind(b,.)
  names(b)[dim(b)[2]]<-" " 
	b_condition<-T
	}
	
	if (n_occ>1 && length(which(gamma==0))!=length(gamma)){
	Gamma<-diag(as.matrix(gamma1))
	StdError<-se[(pfixe+locc+length(Omega)+1):(pfixe+locc+length(Omega)+length(Gamma))] 
	RSE<-cv[(pfixe+locc+length(Omega)+1):(pfixe+locc+length(Omega)+length(Gamma))] 
	b2_condition<-F
	if(length(which(vec==F))!=length(vec)) {
	b2<-data.frame(Gamma,StdError,RSE,row.names=parameters[diag(gamma)!=0]) 
	.<-c(rep("%",(length(Gamma)))) 
	b2<-cbind(b2,.) 
	names(b2)[dim(b2)[2]]<-" "
	b2_condition<-T
	}
	}
	if(n_occ==1 || length(which(gamma==0))!=length(gamma)){b2_condition<-F}
	
	
	
	StdError<-se[(pp-(lSIG-1)):pp]
	RSE<-cv[(pp-(lSIG-1)):pp]
	.<-rep("%",lSIG) 
	#sig.names<-c('sig.inter','sig.slope','sig.interPD','sig.slopePD')[vec2] 
	sig.names<-unlist(c(f[[12]]))[vec2] 
	c<-cbind(data.frame(SIG,StdError,RSE,row.names=sig.names),.)
	names(c)[c(1,dim(c)[2])]<-c("Sigma"," ")
	#boucle sur les doses pour l'affaiche de toues les differents types de reponses
	#e<-list()
	pm<-0
	ll<-lapply(prots,length)
		if(length(which(unlist(lapply(1:nr,function(nr) lapply(1:ll[[nr]], function(ll,prots,nr) length(prots[[nr]][[ll]]),prots=prots,nr=nr)))==0))==0)
    {
	      #e<-lapply(1:nr,function(nr,Dose,subjects,prots) data.frame(subjects,Dose,row.names=as.character(prots[[nr]][1:length(prots[[nr]])])),Dose<-doses[1:length(prots[[nr]])], subjects<-subjects,prots<-prots)
	     e<-lapply(1:nr,function(nr,doses,subjects,prots) data.frame(times=paste(prots[[nr]][1:length(prots[[nr]])]),subjects,doses),Dose<-doses[1:length(prots[[nr]])], subjects<-subjects,prots<-prots)
        pm<-1
        cat("\n","\n")
    }
    
	sink((paste(directory,"\\",output,sep=""))) 
	cat("PFIM 3.2 Option 1",'\n',"\n")
	cat("Project: ", project) 
	cat("\n","\n")
	cat('Date: ', date()) 
	cat("\n","\n")
	cat("\n","\n")
	cat("**************************** INPUT SUMMARY ********************************") 
	cat("\n","\n") 
	cat('Analytical function models : ','\n','\n') 
	ff<-sapply(1:lf, function(lf, form) cat(paste(form[lf]),'\n','\n'), form=form)
	if (b_condition==F) 
  {
  cat('Individual design:','\n') 
  } else
  {
  cat('Population design: ','\n')
  }
 
	for (i in 1:nr){
    if (pm==1){
        cat('Sample times for response:',LETTERS[i],'\n')
        print(e[[i]])
        cat("\n","\n")
    } else{
        cat('Sample times for response:',LETTERS[i],'       Number of subjects per group',' Doses \n')
        cat(paste(prots[[i]],"                                    ",t(t(subjects)),"    ",doses,'\n'),sep="")
        cat("\n","\n")
    }
  }
  
     if (b_condition==T && n_occ>1 && length(which(gamma==0))!=length(gamma)) 
  {
  cat("Number of occasions:",n_occ)
   cat("\n","\n")
  }
  

  cat("Random effect model: Trand = ",Trand) 
  cat("\n","\n")
  
  ff1<-sapply(1:nr,function(nr,sigmainter,sigmaslope) cat(paste('Variance error model response',LETTERS[nr],': (',sigmainter[[nr]],'+',sigmaslope[[nr]],"*f)^2\n")),sigmainter<-sigmainter,sigmaslope<-sigmaslope)
	cat("\n","\n")
  
  if (covariate.model==T || covariate_occ.model==T) {
    cat("Covariate model : ","\n","\n") 
    cat("\t")
    	     if (Trand==1) cat("NB: Covariates are additive on parameters") else cat("NB: Covariates are additive on log parameters")
          cat("\n","\n")
  }
  

	if (covariate.model==T)
	{
        cat("\t")
      	if (covariate_occ.model==T)     cat("Covariates not changing with occasion","\n","\n")
    	     e.cov<-list()
           lcov1<-length(info.covariate[[1]])
           #lcov<-length(info.covariate[[2]][[i]])
          for (i in 1:lcov1){
           #e.cov<-lapply(1:3,function(lcov1,Name,Categories,Proportions,Parameters_associated) data.frame(Name,Categories,Proportions,Parameters_associated),Name<-rep(paste(info.covariate[[1]][[lcov1]]),length(info.covariate[[2]][[lcov1]])), Categories<-info.covariate[[2]][[lcov1]],Proportions<-info.covariate[[3]][[lcov1]],Parameters_associated<-c(paste(info.covariate[[4]][[lcov1]]),rep(".",length(info.covariate[[2]][[lcov1]])-1)))
            Name<-rep(paste(info.covariate[[1]][[i]]),length(info.covariate[[2]][[i]]))
            lcov2<-length(info.covariate[[2]][[i]])
            prem_col<-lapply(1:lcov2,function(lcov2) paste("(",c(lcov2),")",sep=""))
            Categories<-info.covariate[[2]][[i]]
            Proportions<-info.covariate[[3]][[i]]
            References<-c("*",c(rep(" ",length(info.covariate[[2]][[i]])-1)))
            if(length(info.covariate[[4]][[i]])>1) {
                    charc1<-info.covariate[[4]][[i]][1]
                   for (ij in 2: length(info.covariate[[4]][[i]]))
                   {
                      
                      charc<-paste(charc1,info.covariate[[4]][[i]][ij],sep=",")
                      charc1<-charc
                   }
              Parameters_associated<-c(paste(charc1),rep(".",length(info.covariate[[2]][[i]])-1))
            } else 
            {
            Parameters_associated<-c(paste(info.covariate[[4]][[i]]),rep(".",length(info.covariate[[2]][[i]])-1))
            }
    
            e.cov[[i]]<-data.frame(Categories,References,Proportions,row.names=prem_col,check.names=F,check.rows=F)
          cat("\t")
          cat('Covariate',i,':',Name[1],'(',Parameters_associated[1],')',"\n")
                print(e.cov[[i]])
          
          }
          
  }
   
  if (covariate_occ.model==T)
	{  
     
     cat("\n","\n")
     cat("\t") 
     cat("Covariates changing with occasion","\n","\n")
     e.cov1<-list()
     e.cov2<-list()
     
          for (k in 1: length(covariate_occ.name)){
          cat("\t")
          cat ("Covariate ",k,":",covariate_occ.name[[k]],"(",parameter_occ.associated[[k]],")","\n")
            
           lcov1<-length(covariate_occ.category[[k]])
           prem_col1<-lapply(1:lcov1,function(lcov1) paste("(",c(lcov1),")",sep=""))
           Categories<-covariate_occ.category[[k]]
           References<-c("*",c(rep(" ",length(covariate_occ.category[[k]])-1)))
           
           e.cov1[[k]]<-data.frame(Categories,References,row.names=prem_col1,check.names=F,check.rows=F)
           print(e.cov1[[k]])
           cat("\n","\n")
        
          lcov2<-length(sequence.proportions[[k]])
          prem_col2<-lapply(1:lcov2,function(lcov2) paste("(",c(lcov2),")",sep=""))
          Sequences<- list()

           for (sequ in 1:nb_seq[k]){
            Sequences[sequ]<-covariate_occ.sequence[[k]][[sequ]][1]
            for (occ in 2:n_occ){
              Sequences[sequ]<-paste(Sequences[sequ],covariate_occ.sequence[[k]][[sequ]][occ])
            }
           }

          #for (sequ in 1:nb_seq){
          #  Sequences[sequ]<-unlist(lapply(2:n_occ-1, function(n_occ) paste(covariate_occ.sequence[[k]][[sequ]][n_occ],covariate_occ.sequence[[k]][[sequ]][n_occ+1]))) 
          #}
           Proportions<-sequence.proportions[[k]]
          e.cov2[[k]]<-data.frame(unlist(Sequences),Proportions,row.names=prem_col2,check.names=F,check.rows=F)
          #. <- covariate_occ.sequence[[k]]
         
          #e.cov2[[k]]<-cbind(e.cov2[[k]],.,)
          names(e.cov2[[k]])<-c("Sequences","Proportions")
          
          print(e.cov2[[k]])
           cat("\n","\n")
          }
    }
          

	cat("\n","\n")
	cat("Computation of the Fisher information matrix: option = ",option) 
	cat("\n","\n") 
	if (b_condition==F) cat("******************* INDIVIDUAL FISHER INFORMATION MATRIX ******************") else cat("******************* POPULATION FISHER INFORMATION MATRIX ******************") 
	cat("\n","\n") 
	print(unclass(mfisher)) 
	cat("\n","\n") 
	cat("************************** EXPECTED STANDARD ERRORS ************************") 
	cat("\n","\n") 
	cat("------------------------ Fixed Effects Parameters -------------------------") 
  	 cat("\n","\n") 
	print(a) 
	cat("\n","\n") 
	if (b_condition==T){
	cat("------------------------- Variance of Inter-Subject Random Effects ----------------------") 
	cat("\n","\n") 
	print(b) 
	cat("\n","\n") 
  }
  if (n_occ>1 && length(which(gamma==0))!=length(gamma)){
  cat("------------------------- Variance of Inter-Occasion Random Effects ----------------------") 
	cat("\n","\n") 
	print(b2) 
	cat("\n","\n")
	} 
  cat("------------------------ Standard deviation of residual error ------------------------ ")
  cat("\n","\n")
	print(c) 
	cat("\n","\n") 
	cat("******************************* DETERMINANT ********************************") 
	cat("\n","\n") 
	cat(determinant) 
	cat("\n","\n") 
	cat("******************************** CRITERION *********************************") 
	cat("\n","\n") 
	cat(crit) 
	cat("\n","\n")
	summary.ci<-NULL
	summary.power_nni<-NULL
	summary.ci_eq<-NULL
	summary.power_nni_eq<-NULL
	if (covariate.model==T || covariate_occ.model==T)
	{   
			if (covariate.model==F)
	   {    
  
      beta.covariate<-NULL
      npar<-NULL
      }
    	if(compute.power==T || compute.nni==T)
    	{
    	cat("***************************** COMPARISON TEST ********************************") 
    	cat("\n","\n") 
    	summary.ci<-CI.computation(0.05/2,c(unlist(beta.covariate),beta.occ1),se[(p+1):(pfixe+locc)],c(npar,no),Trand)
    	summary.power_nni<-power_nni.computation(alpha/2,given.power,subjects,c(unlist(beta.covariate),beta.occ1),se[(p+1):(pfixe+locc)],c(npar,no),Trand)
    	print(summary.ci)
    	cat("\n","\n")
   	  cat("Type I error =",alpha)
      cat("\n","\n")
      print(summary.power_nni)                                         
    	}
    	cat("\n","\n") 
    	

    	cat("\n","\n") 
    	if(compute.power_eq==T || compute.nni_eq==T)
    	{
    	cat("***************************** EQUIVALENCE TEST ********************************") 
    	cat("\n","\n") 
    	summary.ci_eq<-CI.computation(0.05,c(unlist(beta.covariate),beta.occ1),se[(p+1):(pfixe+locc)],c(npar,no),Trand)
    	summary.power_nni_eq<-power_nni_eq.computation(alpha,given.power,subjects,c(unlist(beta.covariate),beta.occ1),se[(p+1):(pfixe+locc)],c(npar,no),Trand)
    	print(summary.ci_eq)
    	cat("\n","\n")
   	  cat("Type I error =",alpha)
   	  cat("\n")
 	  	cat("Equivalence interval = [log(",paste(exp(interval_eq[1])),"),log(",paste(exp(interval_eq[2])),")]",sep="")
      cat("\n","\n")
      print(summary.power_nni_eq)
      if (length(which(c(unlist(beta.covariate),beta.occ1) < interval_eq[1])) + length(which(c(unlist(beta.covariate),beta.occ1)> interval_eq[2])) > 0)
      cat("\n","NA: Unavailable because beta is not in the equivalence interval (not under H1)")                                          
    	}
    
    	cat("\n","\n") 
    	
	}
	
	    	
	cat("\n","\n")
	d1<-Sys.time()
  print(d1-d) 	
	sink() 
	if(graph.logical==T) {graph(formg_init,prots,graphinf,graphsup,names.datax,names.datay,y.range)} else NULL 
	cat('OK AF cov','\n') 
	if (covariate.model==T)
	   return(list(doses=dose,prots=prots,subjects=subjects,mfisher=mfisher,determinant=determinant,crit=crit,se=se,cv=cv,summary.ci=summary.ci,summary.power_nni=summary.power_nni,summary.ci_eq=summary.ci_eq,summary.power_nni_eq=summary.power_nni_eq))
  if (covariate.model==F)
  	return(list(doses=dose,prots=prots,subjects=subjects,mfisher=mfisher,determinant=determinant,crit=crit,se=se,cv=cv))

}
