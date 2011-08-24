# Author: Martin Fink
# Date: July 2011
#
# Test various possibilities for deriving inverse and traces of symmetric matrices to
# make the PFIM calculation faster and more robust
#
## Input
# none
## Return list
# none
##
rm(list=ls(all=TRUE))

# Program myPFIMv2
names <- dir('../R','model.*.R')
for (iName in names) source(paste('../R/',iName,sep=''))

modelJ <- model.define(model.defaultDEjac())
	a<- model.plot(modelJ,F,F)
	mJac <- model.sens(modelJ)
	mfJac <- model.fish(modelJ, mJac)

## Check out inversion - diagonal of inverted matrix, etc.
mFish <- mfJac[[1]][[1]]
mFinv <- solve(mFish)
sort(eigen(mFish, symm=T, only.values=T)$values)
1/eigen(mFinv, symm=T, only.values=T)$values
# Eigenvalues of inverse are 1/eigen of normal
mEF <- eigen(mFish, symm=T)

# Diagonal elements of A x D x t(A)
#
#			(d1   0    0)		( a11 a21 a31 )
#			(0   d2    0)		( a12 a22 a32 )
#			(0   0    d3)		( a13 a23 a33 )
#
#( a11 a12 a13 )   a11d1  a12d2  a13d3  |  a11a11d1 + a12a12d2 + a13a13d3
#( a21 a22 a23 )   a21d1  a22d2  a23d3  |		a21a21d1 + a22a22d2 + a23a23d3
#( a31 a32 a33 )   a31d1  d32d2  a33d3  |			a31a31d1 + a32a32d2 + a33a33d3
#
diagInvSymm <- function(A) {mEF <- eigen(A, symm=T); return((mEF$vectors*mEF$vectors) %*% (1/mEF$values))}
t(diag(mFinv)/diagInvSymm(mFish))

# Trace of matrix multiplication A x t(B)
#
#			( b11 b21 b31 )
#			( b12 b22 b32 )
#			( b13 b23 b33 )
#( a11 a12 a13 )   a11b11 + a12b12 + a13b13
#( a21 a22 a23 )   	a21b21 + a22b22 + a23b23
#( a31 a32 a33 )   		a31b31 + a32b32 + a33b33
#
sz <- dim(mFish)
B <- matrix(rnorm(sz[1]*sz[2]),sz[1],sz[2]); B <- B + t(B)
#traceMult <- function(A,B) {as.vector(A) %*% as.vector(t(B))} 
traceMult <- function(A,B) {
	n <- length(A)
	dim(A) <- n
	B <- as.matrix(B,byrow=T)
	dim(B) <- n
	return(A %*% B)
} 
system.time(for (i in seq(1e4)) t1 <- traceMult(mFish, B))
system.time(for (i in seq(1e4)) t2 <- sum(diag(mFish %*% B)))
system.time(for (i in seq(1e4)) t3 <- sum(mFish*B))
all.equal(as.numeric(t1),t2)
all.equal(t3,t2)

modelJ <- model.define(model.defaultAE())
Rprof(filename = "output/Rprof.out", append = FALSE, interval = 0.02, memory.profiling=FALSE)
for (i in seq(1e2))
	mfJac <- model.fish(modelJ)
Rprof(NULL); summaryRprof("output/Rprof.out")
