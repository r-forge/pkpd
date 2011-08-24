#PFIM 3.2.2 Multiple responses
#February 2011
#Copyright © PFIM 3.2.2 – Caroline Bazzoli, Thu Thuy Nguyen, Anne Dubois, Sylvie Retout, Emanuelle Comets, France Mentré - Université Paris Diderot – INSERM.
#------------------------------------------------------------------------------------


###########################################################################
# 	Computation of the Population Fisher Information matrix       	 	  #     
#	 							 									 #							     
###########################################################################


#---------------------- PRE - COMPUTATION ------------------------------------
#-----------------------------------------------------------------------------
library(numDeriv)

theta<-beta
#lf<-length(form)
p<-length(theta)

#ld<-length(protPD)                   #meme nombre de protocoles elementaires
#prots<-list(prot,protPD)

prots<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("prot",LETTERS[1:nr],sep=""))
sigmainter<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("sig.inter",LETTERS[1:nr],sep=""))
sigmaslope<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("sig.slope",LETTERS[1:nr],sep=""))
#bornes<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("borne",LETTERS[1:nr],sep=""))

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

lp<-length(prots[[1]])

if (condinit.identical==T) condinit<-rep(condinit[1],lp) else NULL

if (length(which(omega==0))==length(omega)) 
{ls<-length(subjects); subjects<-rep(1,ls)}

#-----------------------------------------------------------------------------
#------------------------ COMPUTE --------------------------------------------


#Function used in Population Fisher information matrix  function 
#------------------------------------
pts<-function(nr,i,prots){
protk<-c()
lprotk<-c()
for (ki in 1:nr){
protk <-c(protk,prots[[ki]][i])  
}
return(protk)
}

#Function used in sensibilite function
#------------------------------------

inter<-function(theta,kit,l)
{
	
	for (i in 1:p) {assign(parameters[i],theta[i]) }
	cond<-eval(condinit[kit])
  res<-lsoda(cond,times=c(time.condinit,theta[p+1]),formED,theta,rtol=RtolEQ,atol=AtolEQ,hmax=Hmax)
  if (nr==1)res<-res[,dim(res)[2]][2]
  else res<-res[,((dim(res)[2]-nr+1):(dim(res)[2]))][,l] [-1]
}


#Sensitivity matrix :
#--------------------
sensibilite<-function(nr,protk,kit,lprotk,condinit)
{	

	for (i in 1:p) {assign(parameters[i],theta[i]) }
	mod<-list()
	dd<-list()
	dd2<-list()
	dtot<-list()
  s<-list()
	s1<-matrix(nrow=p,ncol=0)
  s2dp<-list()
	cond<-eval(condinit[kit])	
	for (l in 1:nr){
	dd2[[l]]<-list()
	dd[[l]]<-list()
	if (lprotk[[l]]==0){
	s<-matrix(nrow=p,ncol=0)
  } else
  {
  mod[[l]]<-lsoda(cond,c(time.condinit,sort(protk[[l]])),formED,theta,rtol=RtolEQ,atol=AtolEQ,hmax=Hmax)
  if (nr==1) mod[[l]]<-mod[[l]][,dim(mod[[l]])[2]][-1]      #single response
  else { mod[[l]]<-mod[[l]][,((dim(mod[[l]])[2]-nr+1):(dim(mod[[l]])[2]))][,l][-1]}
  
  lfd<-lprotk[[l]]
  #computation of gradient and hessian  in one list
  dd[[l]]<-lapply(1:lfd,function(lfd,theta,protk,kit,l) fdHess(c(theta,protk[lfd]),inter,kit=kit,l=l)$gradient,theta=theta,protk=protk[[l]],kit=kit,l=l)
  if (lprotk[[l]]!=0) dd2[[l]]<-lapply(1:lfd,function(lfd,theta,protk,kit,l) hessian(inter,c(theta,protk[lfd]),method="Richardson",kit=kit,l=l),theta=theta,protk=protk[[l]],kit=kit,l=l)
  
  s<-t(matrix(unlist(dd[[l]]),ncol=lprotk[[l]]))[,1:p]
  
  if(lprotk[[l]]==1) s<-matrix(s,nrow=p,ncol=1) else s<-t(s)
  
  }
  if (Trand==2) {
  s<-s*theta
  }
  s1<-cbind(s1,s)
  
  if (lprotk[[l]]!=0) { 
  #matrice hessienne 
  s2d<-NULL
  s2dp[[l]]<-list()
       for (j1 in 1:p)
        {
           for (j in 1:lprotk[[l]])
            {
             if (Trand==2)
             {
             tran<-matrix(dd2[[l]][[j]][,j1][1:p]*theta,ncol=1)
             tran[j1]<-tran[j1]+s[j1,j]/theta[j1]
             s2d<-cbind(s2d,tran)
             }
             else
             s2d<-cbind(s2d,matrix(dd2[[l]][[j]][,j1][1:p],ncol=1))
            
            }
         s2dp[[l]][[j1]]<-s2d
         s2d<-NULL
        }
  }
  }
  #joindre les morceaux de chaque réponse par parametres#
     s3dp<-list()
     for (j4 in 1:p)
     {
     mo1<-NULL
     for (j3 in 1:length(s2dp))
     {
     mo<-NULL
     mo<-s2dp[[j3]][[j4]]
     mo1<-cbind(mo1,mo)
     }
     s3dp[[j4]]<-mo1
     }    
	l1<-list(unlist(mod),s1,s3dp)

	return(l1) 
}


# Fisher matrix for 1 ind in 1 elementary design by block: block A (resA), block B(resB), block C (resC)
#-------------------------------------------------------------------------------------------------------
bloc<-function(protk,sensi,var,invvar,dvar,mod,lprotk,bout4)
{
	resB<-matrix(nrow=p+nr*2,ncol=p+nr*2)     #matrice aléatoire #
	resC<-matrix(rep(0,(p+nr*2)*p),nrow=p+nr*2,ncol=p)     #matrice incomplete option 1 covariance==0#
	resA2<-matrix(ncol=p,nrow=p)
	if (Trand==2) sensia<-sensi/theta 	else sensia<-sensi 
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
	    lprotk<-lapply(1:nr,function(nr,protk) length(protk[[nr]]),protk=protk )
	    lcumprotk<-c(1,cumsum(unlist(lprotk)))
	    if (sum(lcumprotk[-1])>0){
			sensibili<-sensibilite(nr,protk,i,lprotk,condinit)  #i the elementar protocol
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
            for (j6 in 1:length(eltdvar))
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
			if (length(which(dim(var)==0))!=2) {invvar<-solve(as.matrix(var))
      fish<-bloc(protk,sensi,var,invvar,dvar,mod,lprotk,bout4) }
      else {fish<-matrix(c(rep(0,pp*pp)),ncol=pp)}
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
	
return(list(formED=formED,somme=somme,subjects=subjects,prots=prots,se=se,cv=cv,det=det,crit=crit,p=p,pp=pp,xn=xn,sensibili=sensibili)) 

}


graph<-function(prots,graphinf,graphsup,names.datax,names.datay,y.range)
{
		for (i in 1:p) {assign(parameters[i],theta[i]) }
		tt<-list()
		ff1<-list
    e<-list()
    e1<-list()
    tt<-lapply(1:nr, function(nr,graphinf,graphsup) seq(graphinf[[nr]],graphsup[[nr]],0.1),graphinf=graphinf,graphsup=graphsup)
	
		
		if (condinit.identical==T)
		{
				par(mfrow=c(1,nr))
				cond<-eval(condinit[1])

        e<-lapply(1:nr,function(nr,time.condinit,tt) c(time.condinit,tt[[nr]]),time.condinit=time.condinit,tt=tt)

				ff1<-list()
        for (li in 1 : nr){
				  ff1[[li]]<-lsoda(cond,e[[li]],formED,theta,rtol=RtolEQ,atol=AtolEQ,hmax=Hmax)
          if (nr==1) ff1[[li]]<-ff1[[li]][,dim(ff1[[li]])[2]] [-1] else  ff1[[li]]<-ff1[[li]][,((dim(ff1[[li]])[2]-nr+1):(dim(ff1[[li]])[2]))][,li] [-1]
         
        
      
          if (log.logical!=F) plot(tt[[li]],ff1[[li]],type='l',xlab=paste(names.datax[li]),ylab=paste(names.datay[li]),log=log.logical,ylim=c(y.range[[li]],range(ff1[[li]]))[1:2]) else  plot(tt[[li]],ff1[[li]],type='l',xlab=paste(names.datax[li]),ylab=paste(names.datay[li]),ylim=c(y.range[[li]],range(ff1[[li]]))[1:2])
      	
            #ff2[[l]]<-vector(mode="numeric", length(prots[[l]][[i]]))
          ff2<-list()
         for (j in 1:length(prots[[li]]))
				{
				      #ff2[[l]]<-lsoda(cond,c(time.condinit,sort(prots[[l]][[j]])),formED,theta,rtol=RtolEQ,atol=AtolEQ,hmax=Hmax)[,cmpt.of.interest[l]+1][-1]/eval(scal[[l]])
              if (length(prots[[li]][[j]])!=0) {
                ff2[[li]]<-lsoda(cond,c(time.condinit,sort(prots[[li]][[j]])),formED,theta,rtol=RtolEQ,atol=AtolEQ,hmax=Hmax)
                if (nr==1) ff2[[li]]<-ff2[[li]][,dim(ff2[[li]])[2]][-1] else ff2[[li]]<-ff2[[li]][,((dim(ff2[[li]])[2]-nr+1):(dim(ff2[[li]])[2]))][,li][-1]
                 #print(ff2)
				        points(prots[[li]][[j]],ff2[[li]],pch=paste(j),cex=2,col=paste(j))
				         #print(prots[[l]][[j]])
				        
				       }
				       title(paste(names.datay[li],"model ","\n Initial Conditions = ", condinit[1],sep=" "))
         }
			  }     

		}
		else
		{
			par(mfrow=c(nr,length(prots[[1]])))
		  e<-list()
      ff1<-list()
			ff2<-list()	
			for (l in 1:nr){
			 for (ki in 1:length(prots[[l]]))
			 {
        cond<-eval(condinit[ki])
        e<-lapply(1:nr,function(nr,time.condinit,tt) c(time.condinit,tt[[nr]]),time.condinit=time.condinit,tt=tt)
				ff1<-list()
        for (li in 1 : nr){
				  ff1[[li]]<-lsoda(cond,e[[li]],formED,theta,rtol=RtolEQ,atol=AtolEQ,hmax=Hmax)
          if (nr==1) ff1[[li]]<-ff1[[li]][,dim(ff1[[li]])[2]][-1]
          else  ff1[[li]]<-ff1[[li]][,((dim(ff1[[li]])[2]-nr+1):(dim(ff1[[li]])[2]))][,li] [-1]
        }
        if (log.logical!=F) plot(tt[[l]],ff1[[l]],type='l',xlab=paste(names.datax[l]),ylab=paste(names.datay[l]),log=log.logical,ylim=c(y.range[[l]],range(ff1[[l]]))[1:2]) else  plot(tt[[l]],ff1[[l]],type='l',xlab=paste(names.datax[l]),ylab=paste(names.datay[l]),ylim=c(y.range[[l]],range(ff1[[l]]))[1:2])
				
         if (length(prots[[l]][[ki]])!=0) {
         ff2[[l]]<-lsoda(cond,c(time.condinit,sort(prots[[l]][[ki]])),formED,theta,rtol=RtolEQ,atol=AtolEQ,hmax=Hmax)
         if (nr==1) ff2[[l]]<-ff2[[l]][,dim(ff2[[l]])[2]][-1] 
         else  ff2[[l]]<-ff2[[l]][,((dim(ff2[[l]])[2]-nr+1):(dim(ff2[[l]])[2]))][,l][-1]
         points(prots[[l]][[ki]],ff2[[l]],pch=paste(ki),cex=2)
	       
         }
         title(paste(names.datay[l],"model  ","\n Initial Conditions ",  condinit[ki],sep=" "))
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
	pm<-0
	ll<-lapply(prots,length)
	if(length(which(unlist(lapply(1:nr,function(nr) lapply(1:ll[[nr]], function(ll,prots,nr) length(prots[[nr]][[ll]]),prots=prots,nr=nr)))==0))==0){
  	e<-lapply(1:nr,function(nr,subjects,prots) data.frame(times=as.character(prots[[nr]][1:length(prots[[nr]])]),subjects), subjects<-subjects,prots<-prots)
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
	cat('Differential Equations form of the model: ','\n','\n') 
	print(formED)
	cat("\n","\n")
	if (b_condition==T) cat('Population design: ','\n')   else cat('Individual design: ','\n')
  cat("\n","\n")
for (i in 1:nr){
    if (pm==1){
        cat('Sample times for response:',LETTERS[i],'\n')
        print(e[[i]])
    } else{
        cat('Sample times for response:',LETTERS[i],'                         Number of subjects per group \n')
        cat(paste(prots[[i]],"      ", t(t(subjects)),'\n'))
        cat("\n","\n")
    }
  }
  ff1<-sapply(1:nr,function(nr,sigmainter,sigmaslope) cat(paste('Variance error model response',LETTERS[nr],': (',sigmainter[[nr]],'+',sigmaslope[[nr]],"*f)^2\n")),sigmainter<-sigmainter,sigmaslope<-sigmaslope)
	cat("\n","\n")
	cat('Initial Conditions at time', time.condinit,':','\n','\n')
	ff<-sapply(1:lp, function(lp, form) cat(paste(form[[lp]][-1]),'\n'), form=condinit)
	cat("\n","\n") 


  cat("Random effect model: Trand = ",Trand) 
  cat("\n","\n")
  
  ff1<-sapply(1:nr,function(nr,sigmainter,sigmaslope) cat(paste('Variance error model response',LETTERS[nr],': (',sigmainter[[nr]],'+',sigmaslope[[nr]],"*f)^2\n")),sigmainter<-sigmainter,sigmaslope<-sigmaslope)
	cat("\n","\n")
	
		cat("Error tolerance for solving differential equations system: RtolEQ =", RtolEQ,", AtolEQ =", AtolEQ,", Hmax = ",Hmax)
	cat("\n","\n") 

  
	if (b_condition==T) cat("******************* POPULATION FISHER INFORMATION matrix ******************") else  cat("******************* INDIVIDUAL FISHER INFORMATION matrix ******************")
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
	if(graph.logical==T)
  {
  graph(prots,graphinf,graphsup,names.datax,names.datay,y.range) 
  } 
  else NULL 
	cat('OK ODE EVAL',"\n") 	
	return(list(prots=prots,subjects=subjects,mfisher=mfisher,determinant=determinant,crit=crit,se=se,cv=cv))
}



