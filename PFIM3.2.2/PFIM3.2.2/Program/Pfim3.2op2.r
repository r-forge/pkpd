#PFIM 3.2.2 Multiple responses
#February 2011
#Copyright © PFIM 3.2.2 – Caroline Bazzoli, Thu Thuy Nguyen, Anne Dubois, Sylvie Retout, Emanuelle Comets, France Mentré - Université Paris Diderot – INSERM.
#------------------------------------------------------------------------------------


#---------------------- PRE - COMPUTATION ------------------------------------
#-----------------------------------------------------------------------------
library(numDeriv)


theta<-beta
lf<-length(form) #nombre de forme analytique
#nombre de form dans chaque réponse pour savoir quoi prendre en tf
#recover of the model
formg<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("form",LETTERS[1:nr],sep=""))  #sous forme de list
lformg<-lapply(formg,length)  #length(of each model tf or not)
c<-c(0,cumsum(lformg))
doses<-dose
p<-length(theta)
#Dform<-lapply(1:lf,function(p,form,parameters,lf) lapply(1:p,function(p,form,parameters) D(form,parameters[p]), form=form[lf],parameters=parameters),p=p,form=form,parameters=parameters)
#recover the inital data for all the models
prots<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("prot",LETTERS[1:nr],sep=""))
sigmainter<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("sig.inter",LETTERS[1:nr],sep=""))
sigmaslope<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("sig.slope",LETTERS[1:nr],sep=""))
#test pour savoir si IFIM ou PFIM#

#tests si toutes les bornes  ont été specifiée
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
         
         
#test if nom pour les abscisses sont presents pour chaque réponse#  
#if (length(names.datax)!=nr) names.datax<-rep("NULL",nr)
#if (length(names.datay)!=nr) names.datay<-rep("NULL",nr)     
          
#test si graph.inf existe pour chaque réponse#
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
tt1<-lapply(1:nr,function(nr,tt1b) tt1b[[nr]][tt1b[[nr]]>=graphinf[[nr]][1] & tt1b[[nr]]<=graphsup[[nr]][1]],tt1b<-tt1b) 

#vector of dose
if (dose.identical==T) doses<-rep(doses,length(prots[[1]])) else NULL

#evaluation of individual design update number of subjects == 1
b_condition<-F
if (length(which(omega==0))==length(omega)) 
{b_condition<-T; ls<-length(subjects); subjects<-rep(1,ls)}

#-----------------------------------------------------------------------------
#------------------------ COMPUTE --------------------------------------------




#Function used in Population Fisher information matrix  function 
#------------------------------------
pts<-function(nr,i,prots){
protk<-c()
lprotk<-c()
for (k in 1:nr){
protk <-c(protk,prots[[k]][i])  

}
return(protk)
}


#Sensitivity matrix :
#-------------------- 
#Sensitivity matrix :
#--------------------
sensibilite<-function(nr,protk,k,lprotk) #m=modele pk ou pd proti=prelevemet pk et protiPD=prelevementPD#
{	            #protk represente touts les temps de prelevement de toute les reponses pour le k ieme protocle elementaire
  #c<-length(proti) #nombre de prelevement par protocle elementaire
	#c1<-length(protiPD)
	#lprotk<-lapply(1:nr,function(nr,protk,lprotk) length(protk[[nr]]),protk=protk )
  dose<-doses[k]
	s1<-matrix(nrow=p,ncol=0)
	mod<-numeric()
	mod1<-numeric()
  s3d<-list()
	for (i in 1:p) {assign(parameters[i],theta[i]) }
  
	 for (l in 1 :nr){
	    #dat<-se(lprotk,protk,dose,form,tf,l)
	    #mod1<-c(mod1,dat[[1]])
	    #s1<-cbind(s1,dat[[2]])
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
      if (lprotk[[l]]!=0) { 
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
			          if(j1==i) 
                {s2d[[i]][j1,j]<-s2d[[i]][j1,j]+s[i,j]/theta[i]
			          #print(s2d[[i]][j1,j])
			          }
			          }
			          else
			          {
			          s2d[[i]][j1,j]<-eval(D2forma[[vlf[t<=tf[[l]]][1]]][[i]][[j1]])
			          }
			          }  
	         }  
             
	     } 
     
	   mod1<-c(mod1,mod)  
     s1<-cbind(s1,s)
     s3d[[l]]<-s2d
     }
        #joindre les morceaux de chaque réponse par parametres#
    }

     #joindre les morceaux de chaque réponse par parametres#
     s3dp<-list()
    for (j4 in 1:p)
     {
     mo1<-NULL

     for (j3 in 1:length(s3d))
      {
      
     mo<-NULL  
     mo<-s3d[[j3]][[j4]] 
     mo1<-cbind(mo1,mo)
    }
    s3dp[[j4]]<-mo1
     }

     l<-list(mod1,s1,s3dp)
     #print(l)
	return(l)
}


bloc<-function(protk,sensi,var,invvar,dvar,mod,lprotk,bout4)
{
	resB<-matrix(nrow=p+nr*2,ncol=p+nr*2)     #matrice aléatoire #
	resC<-matrix(rep(0,(p+nr*2)*p),nrow=p+nr*2,ncol=p)     #matrice incomplete option 1 covariance==0#
	resA2<-matrix(ncol=p,nrow=p)
	if (Trand==2) sensia<-sensi/theta 	else sensia<-sensi 
	resA1<-(sensia %*% invvar) %*% t(sensia)  #matrice parametres fixe#
	#print(resA1)
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
      #print(resB1)
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
	h1<-cbind(rbind(resA,resC),rbind(t(resC),resB))
	return(h1) 
}

#Population Fisher information matrix
#------------------------------------
fisher<-function()
{	 
	pp<-2*p+2*nr #2 parameters for the residual variance by response#
	somme<-matrix(c(rep(0,pp*pp)),ncol=pp) 
	se<-numeric() 
	cv<-numeric() 
	
	for(i in 1: length(prots[[1]]))  #pour chaque protocole meme nombre de protocole elementaire en pk et pd#
	{
	    protk<-pts(nr,i,prots)
	    print(protk)
	    lprotk<-lapply(1:nr,function(nr,protk) length(protk[[nr]]),protk=protk )
	    
	    lcumprotk<-c(1,cumsum(unlist(lprotk)))
	    if (sum(lcumprotk[-1])>0){
			sensibili<-sensibilite(nr,protk,i,lprotk)  #i the elementar protocol 
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
			
      if(lprotk[[n]]!=0) eltdvar[[n]]<-lapply(1:p, function(p,sensi,t,t1,sigmaslope,n,theta) sigmaslope[[n]]*sensi[p,t:t1]/theta[p],sensi=sensi,t=t,t1=t1,sigmaslope=sigmaslope,n=n,theta=theta)
      }
			t<-NULL
			t1<-NULL
			#element dans dvar  --> diag(c(sigmaslope*(sensi[p,1:lprotk]/theta[p]),sigmaslopePD*(sensi[p,(c+1):(c+c1)]/theta[p])),length(prots))
      #pour reconstituer une liste sur les parametres et non sur le nombre de reponses
      eltdvarf<-list()
      for (j5 in 1:p)
      {
          mo1<-NULL
            for (j6 in 1: length(eltdvar))
            {
              mo<-NULL
              mo<-eltdvar[[j6]][[j5]]
              mo1<-c(mo1,mo)
            }
          eltdvarf[[j5]]<-mo1
      }

      
			dvar<-lapply(1:p,function(dsensi2,bout1,bout4,sigmaslope,sensi,theta,p,protk)
			 			t(dsensi2[[p]])%*%bout1+t(bout1)%*%dsensi2[[p]]+2*diag(unlist(bout4[1:nr]),length(unlist(protk)))*diag(eltdvarf[[p]],length(unlist(protk))),
							dsensi2=dsensi2,bout1=bout1,bout4=bout4,sigmaslope=sigmaslope,sensi=sensi,theta=theta,protk=protk)
  
      var<-t(bout1) %*% sensi + diag(unlist(bout4[1:nr])^2,length(unlist(protk))) 
			#print(bout1) 
      if (length(which(dim(var)==0))!=2) {invvar<-solve(as.matrix(var))
           fish<-bloc(protk,sensi,var,invvar,dvar,mod,lprotk,bout4) }
      else {fish<-matrix(c(rep(0,pp*pp)),ncol=pp)}
      #print(dvar)
      #print(fish)
      }else{fish<-matrix(c(rep(0,pp*pp)),ncol=pp)}
       
			somme<-somme+subjects[i]*fish 	
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
		cv1<-se[1:p]/theta*100 
		cv2<-se[(p+1):(pp-lSIG)]/diag(as.matrix(omega1))*100 
		cv3<-se[(pp-(lSIG-1)):pp]/SIG*100 
		cv<-abs(c(cv1,cv2,cv3))	
		crit<-(det)^(1/pp) 
	}
	
return(list(form=form,somme=somme,subjects=subjects,prots=prots,se=se,cv=cv,det=det,crit=crit,p=p,pp=pp,xn=xn)) 

}




#function to build graph
graph<-function(prots,graphinf,graphsup,names.datax,names.datay,y.range)
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
  		      if(log.logical!=F) plot(tt1[[l]],ff1[[l]],type='l',xlab=paste(names.datax[l]),ylab=paste(names.datay[l]),log=log.logical,ylim=c(y.range[[l]],range(ff1[[l]]))[1:2]) else   plot(tt1[[l]],ff1[[l]],type='l',xlab=paste(names.datax[l]),ylab=paste(names.datay[l]),ylim=c(y.range[[l]],range(ff1[[l]]))[1:2])
        
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
			     for (k in 1:length(prots[[l]]))
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
              ff2[[l]]<-vector(mode="numeric", length(prots[[l]][[k]]))
				      for ( i in 1:length(prots[[l]][[k]]))
				      {
	               dose<-doses[k]
					       t<-prots[[l]][[k]][i]
					       ff2 [[l]][i]<-eval(formg[[l]][t<=tf[[l]]][1])#*dose[k] #				
		          }
           	    points(prots[[l]][[k]],ff2[[l]],pch=paste(k),cex=2)
			         	title(paste("Group",k,"- dose = ", doses[k],sep=" "))
			}
			}
      	
		}		
}

#------------------------------------------------------------------------------
#------------------------------ OUTPUT -----------------------------------

out<-function()
{

  d<-Sys.time()
	f<-fisher() 

	p<-f$p 
	pp<-f$pp 
	se<-f$se 
	cv<-f$cv 
	mfisher<-f$somme 
	determinant<-f$det 
	crit<-f$crit 
	
  StdError<-se[1:p] 
	RSE<-cv[1:p] 
	Beta<-theta 

	a<-data.frame(Beta,StdError,RSE,row.names=parameters) 
	.<-c(rep("%",p)) 
	a<-cbind(a,.) 
  names(a)[dim(a)[2]]<-" "
  
	Omega<-diag(as.matrix(omega1)) 
	StdError<-se[(p+1):(pp-lSIG)] 
	
	RSE<-cv[(p+1):(pp-lSIG)] 
	b_condition<-F
	if(length(which(vec==F))!=length(vec)) {
	b<-data.frame(Omega,StdError,RSE,row.names=parameters[vec]) 
	.<-c(rep("%",(pp-p-lSIG))) 
	b<-cbind(b,.)
  names(b)[dim(b)[2]]<-" "  
	b_condition<-T
	}
	
	StdError<-se[(pp-(lSIG-1)):pp]
	RSE<-cv[(pp-(lSIG-1)):pp]
	.<-rep("%",lSIG) 
	#sig.names<-c('sig.inter','sig.slope','sig.interPD','sig.slopePD')[vec2] 
	sig.names<-unlist(c(f$xn))[vec2] 
	c<-cbind(data.frame(SIG,StdError,RSE,row.names=sig.names),.)
	names(c)[c(1,dim(c)[2])]<-c("Sigma"," ")
	#boucle sur les doses pour l'affaiche de toues les differents types de reponses
	#e<-list()
	pm<-0
	ll<-lapply(prots,length)
	if(length(which(unlist(lapply(1:nr,function(nr) lapply(1:ll[[nr]], function(ll,prots,nr) length(prots[[nr]][[ll]]),prots=prots,nr=nr)))==0))==0){
	   e<-lapply(1:nr,function(nr,Dose,subjects,prots) data.frame(times=as.character(prots[[nr]][1:length(prots[[nr]])]),subjects,Dose),Dose<-doses[1:length(prots[[nr]])], subjects<-subjects,prots<-prots)
	pm<-1
  }
	sink((paste(directory,"\\",output,sep=""))) 
	cat("PFIM 3.2 Option 2",'\n',"\n")
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
  

  cat("Random effect model: Trand = ",Trand) 
  cat("\n","\n")
  
  ff1<-sapply(1:nr,function(nr,sigmainter,sigmaslope) cat(paste('Variance error model response',LETTERS[nr],': (',sigmainter[[nr]],'+',sigmaslope[[nr]],"*f)^2\n")),sigmainter<-sigmainter,sigmaslope<-sigmaslope)
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
	cat("\n","\n")
	d1<-Sys.time()
  print(d1-d) 	
	sink() 
	if(graph.logical==T) {graph(prots,graphinf,graphsup,names.datax,names.datay,y.range)} else NULL 
	cat('OK',"\n") 	
	return(list(doses=doses,prots=prots,subjects=subjects,mfisher=mfisher,determinant=determinant,crit=crit,se=se,cv=cv))
}



