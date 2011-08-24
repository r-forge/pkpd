# Question: Separate respones? Can be easily dealt with as additional times...
# Fisher matrix for 1 ind in 1 elementary design by block: block A (resA), block B(resB), block C (resC)
# avec IOV
#-------------------------------------------------------------------------------------------------------
model.fish.bloc<-function(mod, sensif, obsIdx, model) {
#Prelude
	# Define user-friendly parameters
	nModpar <- model$nModpar
	parModel <- model$parModel
	nObs <- model$nObs

	# mod... model solutions (vector length nTime)
	# sensia .. model sensitivities w.r.t. parameters (matrix nModpar x nTime)
	# obsIdx .. measurement type/observation number for each time point
	
	# Random error variance (additive + proportional)
	sigmainter <- model$parObsErr[1,]; sigmaslope <- model$parObsErr[2,]
	randErr <- 	sigmainter[obsIdx] + sigmaslope[obsIdx]*mod

	# Random effect variance diagonal matrix
	omega <- diag(model$parModelVar)
	
	# Additive (1) or Exponential (2) random effect model?
	# f(x,p,w), P = p+w; P = p*exp(w) => for typical w=0 => df/dP = df/dp
	# df/dw = df/dP * dP/dw; dP/dw = 1 OR w*p*exp(w)
	if (model$parModelVarType %in% c(2,'exp')) { # Multiply sensitivities with parameter
		sensia <- sensif*parModel
	}

	# var = t(omega %*% sensia) %*% sensia + diag(err)
	var<-crossprod(omega%*%sensia, sensia) + diag(randErr^2)

	# As I read again and again... - Don't invert a matrix!!
	#invvar<-try(solve(var))
	#if (is.character(invvar)) {
	#	invvar <- myPseudoinverse(var)
	#	warning("'var' singular. Pseudoinverse used.")
	#}
	mEigen <- eigen(var, symmetric=T)
	invvar <- tcrossprod(mEigen$vectors %*% diag(1/mEigen$values), mEigen$vectors)
	# Try as alternative to use solve(x,y) instead of eigen(x)...*y

#Block A (resA == resA2)
#-------
	resA <- tcrossprod(sensif %*% invvar, sensif)  #matrice parametres fixe#

#Block B
#------- # nOcc = 1
	# resB2 <- all elements of B are traces of matrix multiplications
	#			of symmetric matrices!!
	# sum(diag(ai[[idx]] %*% a1[[jdx]])) = sum(ai[[idx]] * a1[[jdx]])

	resB <- matrix(0,nrow=nModpar+2*nObs, ncol=nModpar+2*nObs)
	a1 <- list(); a2 <- list(); ai <- list();
	# Random error model
	for (idx in seq(nObs)) {
		a1[[idx]]<-2*diag(randErr*(obsIdx==idx))
		a2[[idx]]<-a1[[idx]]*mod
		ai[[idx]]<-(invvar%*%a1[[idx]])%*%invvar
		# Cross-Observation
		for (jdx in 1:idx) {
			resB[nModpar+idx*2-1,nModpar+jdx*2-1]<-sum(ai[[idx]] * a1[[jdx]])
			resB[nModpar+idx*2,  nModpar+jdx*2]  <-sum(diag((a2[[idx]]%*%invvar)%*%(a2[[jdx]]%*%invvar)))
			resB[nModpar+idx*2,  nModpar+jdx*2-1]<-sum(ai[[jdx]] * a2[[idx]])
			if (jdx<idx) resB[nModpar+idx*2-1,nModpar+jdx*2]  <-sum(ai[[idx]] * a2[[jdx]])
		}
	}
	# Variances & random error
	for (i in seq(nModpar)) {
		resB1<-tcrossprod(sensia[i,], sensia[i,])
		resB2<-invvar%*%resB1%*%invvar
		### Is this always the same result?!? Where does index 'i' come into play??
		resB[i,1:i]<-sapply(1:i,function(resB2,sensia,nModpar) sum(resB2%*%sensia[nModpar,] * sensia[nModpar,]),resB2=resB2,sensia=sensia)
		for (idx in 1:nObs) {
			resB[nModpar+idx*2-1,i]<-sum(a1[[idx]] * resB2)
			resB[nModpar+idx*2,i]<-sum(a2[[idx]] * resB2)
		}
	}
	resB<- resB + t(resB)*upper.tri(resB)
	resC<-matrix(0,nrow=nModpar+nObs*2,ncol=nModpar)
	h <- cbind(rbind(resA, resC), rbind(t(resC), resB/2))

	select <- c(parModel, model$parModelVar, as.vector(model$parObsErr)) != 0
	fish <- h[select, select]
	return(list(fish=fish,mod=mod,sensif=sensif,sensia=sensia,var=var))

	# # Random error
	# resB_re <- matrix(nrow=2*nObs, ncol=2*nObs)
	# for (idx in 1:nObs) {
		# a1<-2*diag(randErr[tIndices[idx]:(tIndices[idx+1]-1)])
		# a2<-a1*mod
		# resB_re[idx*2-1,idx*2-1]<-sum(diag(invvar%*%a1%*%invvar%*%a1))
		# resB_re[idx*2,idx*2]<-sum(diag(a2%*%invvar%*%a2%*%invvar))
		# resB_re[idx*2-1,idx*2]<-resB_re[idx*2,idx*2-1]<-sum(diag(ai%*%a2))
	# }
	
	# len<-p
	# resB<-matrix(nrow=p+nr*2,ncol=p+nr*2)     #matrice aléatoire (block B) #
	# resB_int<-resB     #calcul intermédiaire

	# l1<-list()
	# li<-list()
	# l2<-list()
	# l1[[1]]<-a1
	# li[[1]]<-ai
	# l2[[1]]<-a2
	
	
	# if (nr!=1){
		# for (i in 1:(nr-1)){
			# vec<-NULL
			# for (1 in 1:nOcc){
				# vec0<-rep(0,length(unlist(protk)))
				# if(length(protk[[i+1]])!=0){
					# vec0<-c(rep(0,length(protk[[i]])),randErr[[1]][[i+1]])
				# } else {
					# vec0<-vec0
				# }
				# vec<-c(vec,vec0)
			# }
			# a1<-2*diag(vec,length(unlist(protk))*nOcc)
			# ai<-invvar%*%a1%*%invvar
			# resB_int[len+i*2+1,len+i*2+1]<-sum(diag(ai%*%a1))
			# #sigma slope #
			# a2<-a1*mod
			# resB_int[len+i*2+2,len+i*2+2]<-sum(diag(a2%*%invvar%*%a2%*%invvar))
			# resB_int[len+(i*2)+1,len+i*2+2]<-sum(diag(ai%*%a2))
			# resB_int[len+i*2+2,len+(i*2)+1]<-resB_int[len+(i*2)+1,len+i*2+2]
			# l1[[i+1]]<-a1
			# li[[i+1]]<-ai
			# l2[[i+1]]<-a2
		# }
# #sigma inter PD et sigma inter PK#
		# for (i in 1: (nr)){
			# for (j in  2 : (nr)){
				# resB_int[(len+j*2-1),len+i*2-1]<-sum(diag(li[[j]]%*%l1[[i]]))
				# resB_int[len+i*2-1,(len+j*2-1)]<-resB_int[(len+j*2-1),len+i*2-1] #resB[p+1,p+3]<-resB[p+3,p+1]
				# resB_int[len+j*2,len+i*2-1]<-sum(diag(li[[i]]%*%l2[[j]])) 	#resB[p+1,p+4]<-resB[p+4,p+1]
				# resB_int[len+i*2-1,len+j*2]<-resB_int[len+j*2,len+i*2-1]
				
				# resB_int[len+j*2-1,len+i*2]<-sum(diag(li[[j]]%*%l2[[i]]))  #resB[p+2,p+3]<-resB[p+3,p+2]
				# resB_int[len+i*2,len+j*2-1]<-resB_int[len+j*2-1,len+i*2]
				# resB_int[len+j*2,len+i*2]<-sum(diag(invvar%*%l2[[j]]%*%invvar%*%l2[[i]]))
				# resB_int[len+i*2,len+j*2]<-resB_int[len+j*2,len+i*2] #resB[p+2,p+4]<-resB[p+4,p+2]
			# }
		# }
		
		# for (i in 1:len)
		# {
			# resB1<-sensia[i,]%*%t(sensia[i,])
			# resB2<-invvar%*%resB1%*%invvar
			# resB_int[i,1:len]<-sapply(1:len,function(resB2,sensia,len) sum(diag(resB2%*%sensia[len,]%*%t(sensia[len,]))),resB2=resB2,sensia=sensia)
			# resB_int[i,len+1]<-sum(diag(l1[[1]]%*%resB2))
			# resB_int[i,len+2]<-sum(diag(l2[[1]]%*%resB2))
			# for (j in 2:nr) {
				# resB_int[i,len+j*2-1]<-sum(diag(l1[[j]]%*%resB2))
				# resB_int[i,len+j*2]<-sum(diag(l2[[j]]%*%resB2)) 	
			# }
			# resB_int[(len+1):(len+(2*nr)),i]<-resB_int[i,(len+1):(len+(2*nr))]
		# }
	# resB<-resB_int
	# }
	# else {
		# for (i in 1:len)
		# {
			# resB1<-sensia[i,]%*%t(sensia[i,])
			# resB2<-invvar%*%resB1%*%invvar
			# resB_int[i,1:len]<-sapply(1:len,function(resB2,sensia,len) sum(diag(resB2%*%sensia[len,]%*%t(sensia[len,]))),resB2=resB2,sensia=sensia)
			# resB_int[i,len+1]<-sum(diag(l1[[1]]%*%resB2))
			# resB_int[i,len+2]<-sum(diag(l2[[1]]%*%resB2))
			# resB_int[(len+1):(len+2),i]<-resB_int[i,(len+1):(len+2)]	
		# }
		# # Déduire B à partir de B_intermédiare
		# if (nOcc>1 && length(which(gamma==0))!=length(gamma)){
			# resB[1:p,1:p]<-resB_int[1:p,1:p]
			# m<-p+1
			# resB[1:p,m:len]<-nOcc*resB_int[1:p,m:len]
			# resB[m:len,1:p]<-resB[1:p,m:len]
			# resB[1:p,(len+1):(len+2)]<-nOcc*resB_int[1:p,(len+1):(len+2)]
			# resB[(len+1):(len+2),1:p]<-nOcc*resB_int[(len+1):(len+2),1:p]
			# resB[m:len,m:len]<-(nOcc**2)*resB_int[m:len,m:len]
			# resB[m:len,(len+1):(len+2)]<-nOcc*resB_int[m:len,(len+1):(len+2)]
			# resB[(len+1):(len+2),m:len]<-nOcc*resB_int[(len+1):(len+2),m:len]
			# resB[(len+1):(len+2),(len+1):(len+2)]<-resB_int[(len+1):(len+2),(len+1):(len+2)]
		# }
		# else resB<-resB_int
	# }
# # B final
	# resB<-0.5*resB
# # Bloc C
	# resC<-matrix(rep(0,(len+nr*2)*(nModpar+nOcc)),nrow=len+nr*2,ncol=nModpar+nOcc)     #matrice incomplete option 1 covariance==0#
	# # où nModpar est 1 dim de A	

	# select <- c(parModel, model$parModelVar, as.vector(model$parObsErr)) != 0
	# fish <- h[select, select]
	
	return(fish)
}


###############################
#myPseudoinverse <- function (m) 
#{
#	msvd = svd(m)
#	if (length(msvd$d) == 0) {
#		return(array(0, dim(m)[2:1]))
#	}
#	else {
#		return(msvd$v %*% (1/msvd$d * t(msvd$u)))
#	}
#}
