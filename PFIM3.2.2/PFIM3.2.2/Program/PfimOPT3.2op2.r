#PFIM 3.2.2 Multiple responses
#February 2011
#Copyright © PFIM 3.2.2 – Caroline Bazzoli, Thu Thuy Nguyen, Anne Dubois, Sylvie Retout, Emanuelle Comets, France Mentré - Université Paris Diderot – INSERM.
#------------------------------------------------------------------------------------

###########################################################################
# 	Computation of the Population Fisher Information matrix       	  #     
#	 for PK Design with analytical form of the model			 									 #							     
###########################################################################

options(expressions=1000)  
library(numDeriv)			
#---------------------- PRE - COMPUTATION ------------------------------------
#-----------------------------------------------------------------------------
#test pour accepter stdin de la version PFIM 3.0
#test to accept stdin file of version PFIM 3.0
#option object
err2<-tryCatch(get("option"), error=function(e) 4)
if(err2==1) option<-4

#names.datax and names.datay objects
ngraph<-2
vec<-c("x","y")
err1<-tryCatch(names.data_test<-lapply(1:ngraph,function(ngraph,x,vec) get(x), x=paste("names.data",vec[ngraph],sep=""), vec=vec), error=function(e) 4)
 if(length(err1)==4 && err1==4)
 {
 names.datax<-rep("Time",nr)
 names.datay<-rep("Concentration",nr)
 }
#IOV
err4<-tryCatch(get("n_occ"), error=function(e) 4)
if(err4==4) n_occ<-0

#covariate_occ.model
err5<-tryCatch(get("covariate_occ.model"), error=function(e) 4)
if(err5==4) covariate_occ.model<-F

subjects.init<-subjects
subjects.initA<-subjects.init
theta<-beta
lf<-length(form) #nombre de forme analytique
#nombre de form dans chaque réponse pour savoir quoi prendre en tf
#recover of the model
formg<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("form",LETTERS[1:nr],sep=""))  #sous forme de list
lformg<-lapply(formg,length)  #length(of each model tf or not)
c<-c(0,cumsum(lformg))
doses<-dose
p<-length(theta)
omega<-as.matrix(omega)
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
          num[kj]=0
          for (i in which(covariate_occ.sequence[[kj]][[l]]!=covariate_occ.category[[kj]][1])){
            num[kj]=num[kj]+1
            param.asso.occ[[kj]][[l]][[i]]=parameter_occ.associated[[kj]]
            
            for (j  in 1:length(param.asso.occ[[kj]][[l]][[i]])){
              
              occ.effect[[kj]][[l]][[i]][which(param.asso.occ[[kj]][[l]][[i]][j]==parameters)]=1
              beta.occ[[kj]][[l]][[i]][which(param.asso.occ[[kj]][[l]][[i]][j]==parameters)]=beta.covariate_occ[[kj]][[num[kj]]][j]
            } 
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
OMEGA<-diag(c(diag(omega),rep(c(diag(gamma)),pn))) #la matrice OMEGA qui regroupe omega et gamma

#----------- Fin de la partie precomputation pour IOV -------------------------------------------



#Dform<-lapply(1:lf,function(p,form,parameters,lf) lapply(1:p,function(p,form,parameters) D(form,parameters[p]), form=form[lf],parameters=parameters),p=p,form=form,parameters=parameters)
#recover the inital data for all the models
if(identical.times==T){
prots<-lapply(1:nr,function(nr,x) get(x),x=paste("prot",LETTERS[1],sep=""))
}else
prots<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("prot",LETTERS[1:nr],sep=""))
 
Q<-length(prots[[1]])
#l_prote<-lapply(1:nr,function(p,forma,parameters,lf) lapply(1:p,function(p,forma,parameters) D(forma,parameters[p]), forma=forma[lf],parameters=parameters),p=p,forma=forma,parameters=parameters)
#nombre de prelevement par prot elementaire
#l_prote<-lapply(1:nr,function(nr,Q,prots) lapply(1:Q, function(nr,prots,Q) length(prots[[Q]]),Q=Q,prots=prots[[nr]]),Q=Q,prots=prots)
#somme des prevement spa rpto selon les reponses
l_prote<-list()
l2<-list()
for (l in 1:Q){
  l_prote[[l]]<-0
 
  for (i in 1:nr){
   l_prote[[l]]<-sum(l_prote[[l]],length(prots[[i]][[l]]))
 #l2[[l]]<-sum(l2[[l]],l_prote[[i]][[l]])
  
   l2[[i]]<-lapply(1:Q, function(Q,prots) length(prots[[Q]]),prots=prots[[i]])
   
 }
}

#lapply(1:Q,function(nr,l_prote,Q) lapply(1:nr, function(Q,l_prote,nr) sum(l_prote[[nr]][[Q]]),nr=nr,l_prote=l_prote),nr=nr,l_prote=l_prote)

sigmainter<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("sig.inter",LETTERS[1:nr],sep=""))
sigmaslope<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("sig.slope",LETTERS[1:nr],sep=""))

#tests si toutes les bornes
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
          
                    


#algo simplex
if (algo.option=="SIMP")  {
if (identical.times==T){
upper<-lapply(1:nr,function(nr,x) get(x),x=paste("upper",LETTERS[1],sep=""))
lower<-lapply(1:nr,function(nr,x) get(x),x=paste("lower",LETTERS[1],sep="")) }
else{
upper<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("upper",LETTERS[1:nr],sep=""))
lower<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("lower",LETTERS[1:nr],sep=""))
}
}

if (algo.option=="FW")  {
if (identical.times==T){
nwind<-lapply(1:nr,function(nr,x) get(x),x=paste("nwind",LETTERS[1],sep=""))
sampwin<-lapply(1:nr,function(nr,x) get(x),x=paste("sampwin",LETTERS[1],sep=""))
#sampwin<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("sampwin",LETTERS[1:nr],sep=""))
nsamp<-lapply(1:nr,function(nr,x) get(x),x=paste("nsamp",LETTERS[1],sep=""))
nmaxpts<-lapply(1:nr,function(nr,x) get(x),x=paste("nmaxpts",LETTERS[1],sep=""))
nminpts<-lapply(1:nr,function(nr,x) get(x),x=paste("nminpts",LETTERS[1],sep=""))
}else{
nwind<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("nwind",LETTERS[1:nr],sep=""))
sampwin<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("sampwin",LETTERS[1:nr],sep=""))
#sampwin<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("sampwin",LETTERS[1:nr],sep=""))
nsamp<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("nsamp",LETTERS[1:nr],sep=""))
nmaxpts<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("nmaxpts",LETTERS[1:nr],sep=""))
nminpts<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("nminpts",LETTERS[1:nr],sep=""))
}
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
  }     
} 


#attrape si il  ya une erreur#
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

#allow to build graph 
tt1<-list()
		#tt<-lapply(1:nr,function(nr,graphinf,graphsup) seq(graphinf,graphsup,0.1),graphinf<-graphinf[[nr]],graphsup<-graphsup[[nr]])	
	  for (l in 1:nr){
          tt1[[l]]<-list()
          if (length(bornes[[l]])==1){
            #if( bornes[[l]][[1]][1]> graphsup[[l]]) stop(" The upper limit for the graph  ",l,"  is too small")
            tt1[[l]][[1]]<-c(bornes[[l]][[1]][1]:graphsup[[l]])}
          else {
  		      for (i in 1: (length(bornes[[l]])-1)){
  		        tt1[[l]][[i]]<-c(bornes[[l]][[i]][1]:bornes[[l]][[i]][2])
  		        if(i==length(bornes[[l]])-1){
  		          #if (bornes[[l]][[i]][2]>graphsup[[l]]) stop(" The upper limit (graph.sup)  for the response  ",LETTERS[ l]  ,"  is too small")
  		          tt1[[l]][[i+1]]<-c(bornes[[l]][[i]][2]:graphsup[[l]])}
             }
          }
    }
tt1b<-lapply(1:nr,function(nr,tt1) unlist(tt1[[nr]]),tt1<-tt1)
tt1<-lapply(1:nr,function(nr,tt1b) tt1b[[nr]][tt1b[[nr]]>=graphinf[[nr]][1] & tt1b[[nr]]<=graphsup[[nr]][1]],tt1b<-tt1b) 

#test if we have NULL in the prots in this case we can not optimise with same sampling time 
ll<-lapply(prots,length)
		if(identical.times==T && length(which(unlist(lapply(1:nr,function(nr) lapply(1:ll[[nr]], function(ll,prots,nr) length(prots[[nr]][[ll]]),prots=prots,nr=nr)))==0))!=0)
		       identical.times=F

if (subjects.input==1) 
	{
	Nt<-lapply(1:nr, function(nr,subjects.init,l2) sum(subjects.init*unlist(l2[[nr]])),subjects.init=subjects.init,l2=l2)
	Ntot<-sum(unlist(Nt))
	subjects.total<-sum(subjects.init)
	subjects.initN<-subjects.init
	subjects.init<-subjects.init/subjects.total
	subjects<-subjects.init
	} else    
	
subjects<-subjects.init


if (dose.identical==T) doses<-rep(doses,length(prots[[1]])) else NULL
condinit<-NULL

#-----------------------------------------------------------------------------
#------------------------ COMPUTE --------------------------------------------

#Function used in Population Fisher information matrix  function 
#------------------------------------
pts<-function(nr,i,prots){
protk<-c()
lprotk<-c()
for (k in 1:nr){
protk <-c(protk,prots[[k]][i])  
#lprotk<-c(protk,prots[[k]][i])
}
return(protk)
}



#Sensitivity matrix :
#--------------------
#sensibilite<-
#nr,protk,pfixe,p,kjk,lprotk,doses,categories.cov,names.cov,beta.covariate,theta,tab.cov.param
sensibilite<-function(nr,protk,kjk,lprotk,doses) #m=modele pk ou pd proti=prelevemet pk et protiPD=prelevementPD#
{	            #protk represente touts les temps de prelevement de toute les reponses pour le k ieme protocle elementaire
  #c<-length(proti) #nombre de prelevement par protocle elementaire
	#c1<-length(protiPD)
	#lprotk<-lapply(1:nr,function(nr,protk,lprotk) length(protk[[nr]]),protk=protk )
	dose<-doses[kjk]
	s1<-matrix(nrow=p,ncol=0)
	mod<-numeric()
	mod1<-numeric()
  s3d<-list()
	for (i in 1:p) {assign(parameters[i],theta[i]) }

	 for (l in 1 :nr){

	    mod<-numeric()
      s<-matrix(nrow=p,ncol=lprotk[[l]])    
      s2d<-list() 
      s4d<-list()
      vlf<-1:lformg[[l]]
      forma<-formg[[l]]
      lf<-length(forma)
      #derivées premieres par reponses
      Dforma<-lapply(1:lf,function(p,forma,parameters,lf) lapply(1:p,function(p,forma,parameters) D(forma,parameters[p]), forma=forma[lf],parameters=parameters),p=p,forma=forma,parameters=parameters)
      #derivées secondes par reponses et par parametres
      lf2<-length(unlist(Dforma))
      p1<-p
      D2forma<-lapply(1:lf,function(p,Dforma,parameters,lf) lapply(1:p,function(p,Dforma,parameters) lapply(1:p1, function(p1,Dforma,parameters) D(Dforma,parameters[p1]),Dforma=Dforma[[p]],parameters=parameters),Dforma=Dforma[[lf]],parameters=parameters),p=p,Dforma=Dforma,parameters=parameters)
			#print(D2forma)    
      for (i in 1:p){
      s2d[[i]]<-matrix(nrow=p,ncol=lprotk[[l]])
	    s4d[[i]]<-matrix(nrow=p,ncol=lprotk[[l]]) 
             
	         for (j in 1:lprotk[[l]])
		       {   
               t<-protk[[l]][j]
			         ff1<-formg[[l]][t<=tf[[l]]][1]	     
			         u<-eval(Dforma[[vlf[t<=tf[[l]]][1]]][[i]])#*dose # 
               if (Trand ==2)  u<-theta[i]*u
			         mod[j]<-eval(ff1)#*dose #
               s[i,j]<-u
               
               for ( j1 in 1:p)
               {
			          if (Trand==2)
			          {
			          int<-theta[j1]
			          s4d[[i]][j1,j]<-eval(D2forma[[vlf[t<=tf[[l]]][1]]][[i]][[j1]])
			          s2d[[i]][j1,j]<-eval(D2forma[[vlf[t<=tf[[l]]][1]]][[i]][[j1]])*int#*dose#
         		
			          if(j1==i) s2d[[i]][j1,j]<-s2d[[i]][j1,j]+s[i,j]/theta[i]
			          }
			          else
			          {
			          s2d[[i]][j1,j]<-eval(D2forma[[vlf[t<=tf[[l]]][1]]][[i]][[j1]])
			          }
			          }  
	         }  
              
	     } 
     #print(s4d) 
	   mod1<-c(mod1,mod)  
     s1<-cbind(s1,s)
     s3d[[l]]<-s2d
     }
     #joindre les morceaux de chaque réponse par parametres#
     s3dp<-list()
     for (j4 in 1:p)
     {
     mo1<-NULL
     for (j3 in 1:nr)
     {
     mo<-NULL
     mo<-s3d[[j3]][[j4]]
     mo1<-cbind(mo1,mo)
     }
     s3dp[[j4]]<-mo1
     }
     l<-list(mod1=mod1,s1=s1,s3dp=s3dp)
     #print(l)
	   
	return(l)
}


# Fisher matrix for 1 ind in 1 elementary design by block: block A (resA), block B(resB), block C (resC)
#-------------------------------------------------------------------------------------------------------
#bloc<-function(protk,sensi,var,invvar,dvar,mod,lprotk,bout4)
#protk,sensi,sensi.rand,var,invvar,pfixe,p,dvar,mod,lprotk,bout4,n.theta,categories.cov,sensi.theta
bloc<-function(protk,sensi,var,invvar,dvar,mod,lprotk,bout4)
{
	resB<-matrix(nrow=p+nr*2,ncol=p+nr*2)     #matrice aléatoire #
	resC<-matrix(rep(0,(p+nr*2)*p),nrow=p+nr*2,ncol=p)     #matrice incomplete option 1 covariance==0#
	resA2<-matrix(ncol=p,nrow=p)
	if (Trand==2) sensia<-sensi/theta 	else sensia<-sensi #?????#
	resA1<-(sensia %*% invvar) %*% t(sensia)  #matrice parametres fixe#
	
	l1<-list()
	li<-list()
	l2<-list()
	#premiere réponse
	c<-cumsum(unlist(lprotk))
	long<-c[nr]-c[1]
  a1<-2*diag(c(bout4[[1]],rep(0,long)),length(unlist(protk)))
	ai<-invvar%*%a1%*%invvar
	resB[p+1,p+1]<-sum(diag(ai%*%a1))
	a2<-a1*mod
	resB[p+2,p+2]<-sum(diag(a2%*%invvar%*%a2%*%invvar))
	resB[p+1,p+2]<-sum(diag(ai%*%a2))
	resB[p+2,p+1]<-resB[p+1,p+2]
	l1[[1]]<-a1
	li[[1]]<-ai
	l2[[1]]<-a2
	
	if (nr!=1){
	  for (i in 1:(nr-1))
   {
	    long<-c[i]
      vec<-rep(0,length(unlist(protk)))
      if(length(protk[[i+1]])!=0){
        vec[long+1:(lprotk[[i+1]])]=bout4[[i+1]]
      }
      else{
        vec<-vec
      } 
  
	   a1<-2*diag(vec,length(unlist(protk)))
	   ai<-invvar%*%a1%*%invvar
	   resB[p+i*2+1,p+i*2+1]<-sum(diag(ai%*%a1))
	   a2<-a1*mod
      resB[p+i*2+2,p+i*2+2]<-sum(diag(a2%*%invvar%*%a2%*%invvar))
      resB[p+(i*2)+1,p+i*2+2]<-sum(diag(ai%*%a2))
      resB[p+i*2+2,p+(i*2)+1]<-resB[p+(i*2)+1,p+i*2+2]
      l1[[i+1]]<-a1
	     li[[i+1]]<-ai
	     l2[[i+1]]<-a2
  }
 
  for (i in 1: (nr)){
    for (j in  2 : (nr)){

	     resB[(p+j*2-1),p+i*2-1]<-sum(diag(li[[j]]%*%l1[[i]])) 
       resB[p+i*2-1,(p+j*2-1)]<-resB[(p+j*2-1),p+i*2-1] #resB[p+1,p+3]<-resB[p+3,p+1]
	     resB[p+j*2,p+i*2-1]<-sum(diag(li[[i]]%*%l2[[j]])) 	#resB[p+1,p+4]<-resB[p+4,p+1]
       resB[p+i*2-1,p+j*2]<-resB[p+j*2,p+i*2-1]
	     resB[p+j*2-1,p+i*2]<-sum(diag(li[[j]]%*%l2[[i]]))  #resB[p+2,p+3]<-resB[p+3,p+2]
	     resB[p+i*2,p+j*2-1]<-resB[p+j*2-1,p+i*2]
	     resB[p+j*2,p+i*2]<-sum(diag(invvar%*%l2[[j]]%*%invvar%*%l2[[i]]))
       resB[p+i*2,p+j*2]<-resB[p+j*2,p+i*2] #resB[p+2,p+4]<-resB[p+4,p+2]
	}
	}
	
	for (i in 1:p)
	{
		resA3<-invvar%*%dvar[[i]]%*%invvar
		resA2[i,]<-sapply(1:p,function(dvar,resA3,p) sum(diag(resA3%*%dvar[[p]])),dvar=dvar,resA3=resA3)
    resB1<-sensi[i,]%*%t(sensi[i,])
		resB2<-invvar%*%resB1%*%invvar
		resB[i,1:p]<-sapply(1:p,function(resB2,sensi,p) sum(diag(resB2%*%sensi[p,]%*%t(sensi[p,]))),resB2=resB2,sensi=sensi)
		resB[i,p+1]<-sum(diag(l1[[1]]%*%resB2))
		resB[i,p+2]<-sum(diag(l2[[1]]%*%resB2))

    for ( j in 2:nr) 
    {
		  resB[i,p+j*2-1]<-sum(diag(l1[[j]]%*%resB2))
		  resB[i,p+j*2]<-sum(diag(l2[[j]]%*%resB2)) 	
    }
      resB[(p+1):(p+(2*nr)),i]<-resB[i,(p+1):(p+(2*nr))]
  
  #BLOCKc
    for( j in 1:p)
    {
      resC1<-invvar%*%dvar[[j]]%*%invvar
		  resC[i,j]<-sum(diag(resC1%*%resB1))
		  resC[p+1,j]<-sum(diag(resC1%*%l1[[1]]))
		  resC[p+2,j]<-sum(diag(resC1%*%l2[[1]]))
		  for ( j2 in 2:nr) {  
		    resC[p+j2*2-1,j]<-sum(diag(resC1%*%l1[[j2]]))
		    resC[p+j2*2,j]<-sum(diag(resC1%*%l2[[j2]]))
      }
    }
  
  }
  } 
  else
  {
        #cas une seule reponse
        for (i in 1:p)
      	{
   	      resA3<-invvar%*%dvar[[i]]%*%invvar
          resA2[i,]<-sapply(1:p,function(dvar,resA3,p) sum(diag(resA3%*%dvar[[p]])),dvar=dvar,resA3=resA3)
          resB1<-sensi[i,]%*%t(sensi[i,])
          resB2<-invvar%*%resB1%*%invvar
		      resB[i,1:p]<-sapply(1:p,function(resB2,sensi,p) sum(diag(resB2%*%sensi[p,]%*%t(sensi[p,]))),resB2=resB2,sensi=sensi)
		      resB[i,p+1]<-sum(diag(l1[[1]]%*%resB2))
		      resB[i,p+2]<-sum(diag(l2[[1]]%*%resB2))
		      resB[(p+1):(p+2),i]<-resB[i,(p+1):(p+2)]	
	       #matrice bloc C
	         for(j in 1:p)
          {
            resC1<-invvar%*%dvar[[j]]%*%invvar
		        resC[i,j]<-sum(diag(resC1%*%resB1))
		        resC[p+1,j]<-sum(diag(resC1%*%l1[[1]]))
		        resC[p+2,j]<-sum(diag(resC1%*%l2[[1]]))
		       
		      }
		    }
  }
  
 	resA<-resA1+0.5*resA2
	resB<-0.5*resB
	resC<-0.5*resC
	h<-cbind(rbind(resA,resC),rbind(t(resC),resB))

	return(h) 
}

#prots,subjects,condinit,doses,Ntot,ind=0,l_prote,nq,categories,names.cov,proportions.cov,beta.covariate,tab_mod_prop,pfixe,p,theta,tab.cov.param
fisher<-function(prots,subjects,condinit,doses,Ntot,ind=0,l_prote,nq,categories,names.cov,proportions.cov,beta.covariate,tab_mod_prop,pfixe,p,theta2,tab.cov.param,n_occ,lposs)
{	
  
  #Nbsubjects<-Ntot/sum((nq+np)*subjects) 
  Nbsubjects<-Ntot/sum((unlist(l_prote)*subjects)) 
  #n<-np+nq 
  n<-sum(unlist(l_prote))              #nq=nombre de prelevements par protocoles elementaires#
	if(ind==2) Nbsubjects<-1
	pp<-2*p+2*nr #2 parameters for the residual variance by response#
	somme<-matrix(c(rep(0,pp*pp)),ncol=pp) 
	se<-numeric() 
	cv<-numeric() 
	
	#for(i in 1: length(prots[[1]])) 
  for(i in 1: length(l_prote))  #pour chaque protocole meme nombre de protocole elementaire en pk et pd#
	{

      protk<-pts(nr,i,prots)
	    lprotk<-lapply(1:nr,function(nr,protk) length(protk[[nr]]),protk=protk )
	    lcumprotk<-c(1,cumsum(unlist(lprotk)))
	    if (lcumprotk[2]>0){
			sensibili<-sensibilite(nr,protk,i,lprotk,doses)  #i the elementar protocol
			mod<-sensibili[[1]] 
			sensi<-sensibili[[2]] 
      dsensi2<-sensibili[[3]]
			bout1<-omega%*%sensi
			bout4<-list()
			t1<-0
			t<-0
			eltdvar<-list()
			for (n in 1:(nr)){
			t<-1+t1
			t1<-t-1+lprotk[[n]]
			if(lprotk[[n]]==0){
			bout4[[n]]<-NULL 
      if(n==nr){bout4[[n+1]]=Inf}
      } else {
			bout3<-sigmainter[[n]]+sigmaslope[[n]]*mod[t:t1]
			bout4[[n]]<-bout3
			}
			eltdvar[[n]]<-lapply(1:p, function(p,sensi,t,t1,sigmaslope,n,theta) sigmaslope[[n]]*sensi[p,t:t1]/theta[p],sensi=sensi,t=t,t1=t1,sigmaslope=sigmaslope,n=n,theta=theta)
      }
			t<-NULL
			t1<-NULL
			#element dans dvar  --> diag(c(sigmaslope*(sensi[p,1:lprotk]/theta[p]),sigmaslopePD*(sensi[p,(c+1):(c+c1)]/theta[p])),length(prots))
      #pour reconstituer une liste sur les parametres et non sur le nombre de reponses
      eltdvarf<-list()
      for (j5 in 1:p)
      {
          mo1<-NULL
            for (j6 in 1:nr)
            {
              mo<-NULL
              mo<-eltdvar[[j6]][[j5]]
              mo1<-c(mo1,mo)
            }
          eltdvarf[[j5]]<-mo1
      }
			dvar<-lapply(1:p,function(dsensi2,bout1,bout4,sigmaslope,sensi,theta,p,protk)
			 			t(dsensi2[[p]])%*%bout1+t(bout1)%*%dsensi2[[p]]+2*diag(unlist(bout4),length(unlist(protk)))*diag(eltdvarf[[p]],length(unlist(protk))),
							dsensi2=dsensi2,bout1=bout1,bout4=bout4,sigmaslope=sigmaslope,sensi=sensi,theta=theta,protk=protk)
      var<-t(bout1) %*% sensi + diag(unlist(bout4[1:nr])^2,length(unlist(protk))) 

			if (length(which(dim(var)==0))!=2) {invvar<-solve(as.matrix(var))
           fish<-bloc(protk,sensi,var,invvar,dvar,mod,lprotk,bout4) }
      else {fish<-matrix(c(rep(0,pp*pp)),ncol=pp)} 
      }else{fish<-matrix(c(rep(0,pp*pp)),ncol=pp)}
      
			somme<-somme+ subjects[i]*Nbsubjects*fish 

		}
	
	
	#SIG<<-c(sig.inter,sig.slope,sig.interPD,sig.slopePD)
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
	
	omega<-as.matrix(omega)
	vec<<-diag(omega)!=0 
	omega1<<-as.matrix(as.matrix(omega[,vec])[vec,]) 
	vec<-c(rep(T,p),vec,vec2) 
 		#remove the row i and the col i of somme if omega[i,i]=0  
	somme<-somme[,vec][vec,] 
	if (ind==2) l<-list(somme) else
    {
	   pp<-dim(somme)[1] 
	   det<-det(somme)	 #Fisher determinant
#	   det<-crit(criterion,somme,pp)
	   if (det==0) crit<-1e-99	else crit<-(det)^(1/pp)	  #dans le cas où le protocole est degenere
     invcrit<-1/crit 
	  
     if (ind==0) l<-invcrit else
		    {
		    f<-fisher.cv(somme,pp,SIG,omega1)
		    inv<-f[[1]]
		    se<-f[[2]]
		    cv<-f[[3]]
		    l<-list(somme=somme,doses=dose,prots=prots,inv=inv,se=se,cv=cv,det=det,crit=crit,p=p,pp=pp,subjects=subjects,xn=xn)
 		     }
 	}
 	
 return(l)
  
}

fisher.cv<-function(somme,pp,SIG,omega1)
{
	inv<-try(solve(as.matrix(somme)))
	if(!is.null(attributes(inv)$class))
	{
		se<-rep(NA,pp)
		cv<-se
	}		
	else
		{
		se<-sqrt(diag(inv)) 
		cv1<-se[1:p]/theta*100 
		cv2<-se[(p+1):(pp-lSIG)]/diag(as.matrix(omega1))*100 
		cv3<-se[(pp-(lSIG-1)):pp]/SIG*100 
		cv<-abs(c(cv1,cv2,cv3))
		}
	l<-list(inv,se,cv)
	return(l)
}

#function to build graph
graph<-function(prots,protopti,graphinf,graphsup,names.datax,names.datay,y.range)
{
		for (i in 1:p) {assign(parameters[i],theta[i]) }

		if ((dose.identical==T)||(dose.identical==F && algo.option=="FW"))
		{
		    dose<-doses[1] 
		    ff1<-list()
		    ff2<-list()
        par(mfrow=c(1,nr))
        for (l in 1 : nr){
		        ff1[[l]]<-vector(mode="numeric", length(tt1[[l]]))
				    for (j in 1:length(tt1[[l]]))
	          {
				      t<-tt1[[l]][j]
				      ff1[[l]][j]<-eval(formg[[l]][t<=tf[[l]]][1])
	
  		      }
  		      if(log.logical!=F) plot(tt1[[l]],ff1[[l]],type='l',xlab=paste(names.datax[l]),ylab=paste(names.datay[l]),log=log.logical,ylim=c(y.range[[l]],range(ff1[[l]]))[1:2]) else   plot(tt1[[l]],ff1[[l]],type='l',xlab=names.datax[l],ylab=paste(names.datay[l]),ylim=c(y.range[[l]],range(ff1[[l]]))[1:2])
        
            #points for the design on the graph
				          
				          for ( i in 1:length(protopti[[l]]))
				          {
				          ff2[[l]]<-vector(mode="numeric", length(protopti[[l]][[i]]))
				            for ( j in 1:length(protopti[[l]][[i]])) {
				    	          t<-protopti[[l]][[i]][j]
					              ff2[[l]][j]<-eval(formg[[l]][t<=tf[[l]]][1])#*dose[1] #
		                }
				                #points(prots[[l]][[i]],ff2[[l]],pch=paste(i),cex=2,col=paste(i))
				                points(protopti[[l]][[i]],ff2[[l]],pch=paste(i),cex=2,col=paste(i))
	                }
             
		      }
      }
		else
		{
	    ff1<-list()
	    par(mfrow=c(nr,length(protopti[[1]])))
			for (l in 1 : nr){
			  
			 ff1[[l]]<-vector(mode="numeric", length(tt1[[l]]))
			     for (k in 1:length(protopti[[l]]))
			     {
			#windows()
			#par(mfrow=c(1,2))
			       dose<-doses[k]
				      for (j in 1:length(tt1[[l]]))
				      {
				      	 t<-tt1[[l]][j]
				       ff1[[l]][j]<-eval(formg[[l]][t<=tf[[l]]][1])#*dose[k]  #
                
			      }
	         
			      if (log.logical!=F) 
			          plot(tt1[[l]],ff1[[l]],type='l',xlab=paste(names.datax[l]),ylab=paste(names.datay[l]),log=log.logical,ylim=c(y.range[[l]],range(ff1[[l]]))[1:2]) else plot(tt1[[l]],ff1[[l]],type='l',xlab=paste(names.datax[l]),ylab=paste(names.datay[l]),ylim=c(y.range[[l]],range(ff1[[l]]))[1:2])
				
        
              ff2<-list()
              ff2[[l]]<-vector(mode="numeric", length(protopti[[l]][[k]]))
				      for ( i in 1:length(protopti[[l]][[k]]))
				      {
	               dose<-doses[k]
					       t<-protopti[[l]][[k]][i]
					       ff2 [[l]][i]<-eval(formg[[l]][t<=tf[[l]]][1])#*dose[k] #				
		          }
           	    #points(prots[[l]][[k]],ff2[[l]],pch=paste(k),cex=2,col=paste(k))
           	    points(protopti[[l]][[k]],ff2[[l]],pch=paste(k),cex=2,col=paste(k))
			         	title(paste("Group",k,"- dose = ", doses[k],sep=" "))
			}
			}
      	
		}		
}






