##' Evaluate the sensitivities (gradients) of a PFIM model
##'
##' Evaluate the sensitivities of a PFIM model object with respect to
##' the model parameters.
##' @title Evaluate sensitivities (gradients)
##' @param model a PFIM model object
##' @param time an optional numeric vector of times at which to evaluate
##' @param useHess logical - whether to use the \code{fdHess} function from
##'    the \pkg{nlme} package or the \code{gradient} function from the
##'    \pkg{numDeriv} package.
##' @param \dots optional arguments, none used at present
##' @return a three-dimensional array of model function values and sensitivities
##' 
sensitivity <-
    function(model, time=model$time, useHess=FALSE, ...)
    UseMethod("sensitivity")

sensitivity.PFIMmodAE <-
    function(model, time=model$time, useHess=FALSE, ...) {
    if (is.numeric(res <-
                   try(solve(model, time, calcSens=TRUE),
                       silent=TRUE)))
        return(res)
    sensNum(model, time, useHess)
}

sensitivity.PFIMmodDE <-
    function(model, time=model$time, useHess=FALSE, ...) {
    if (!is.null(model$JacIC)) return(sensJac(model, time))
    sensNum(model, time, useHess)
}

sensNum <- function(model, time=model$time, useHess=FALSE) {
    .model.solve.jac.interface <- function(parModel, model, time) {
        model$parModel <- parModel
        res <- solve(model, time)
        ## Return a vector for calculating Jacobian
        as.vector(res)
    }

    .model.solve.fdHess.interface <- function(parModel, model, t, iObs, iArms) {
        model$parModel <- parModel
        res <- solve(model, c(model$tRange[1], t))
        res[2,iObs+(iArms-1)*model$nObs,1]
    }

                                        # Define user-friendly parameters
    nModpar <- model$nModpar
    nObs <- model$nObs
    nArms <- model$nArms
    parModel <- model$parModel
    parArms <- model$parArms
    
                                        # Options for jacobian
    methods.args=list(r=4,eps=1e-4,d=1e-2,zero.tol=sqrt(.Machine$double.eps))
    ## If higher accuracy of Jacobian is wanted, increase from r=4 to r=6, eps=1e-4
    
    res <- array(dim=c(length(time), nObs*nArms, 1+nModpar))
    res[,,1] <- solve(model, time)
    
    ## Two different ways to derive the partial derivatives:
                                        # jacobian - library(numDeriv)
                                        # fdHess - library(nlme) - used by PFIM
    if (!useHess) {                     # Jacobian
        require(numDeriv)
        df1<-jacobian(.model.solve.jac.interface,
                      model$parModel, model=model, time=time,
                      method.args=methods.args)
        for (iModpar in seq(nModpar)) {
            tmp <- df1[,iModpar]
            dim(tmp) <- c(length(time),model$nObs*nArms)
            res[,,1+iModpar] <- tmp
        }
    } else { # fdHess
        require(nlme)
        for (it in seq(time)) {
            for (iObs in seq(nObs)) {
                for (iArms in seq(nArms)) {
                    res[it,iObs+(iArms-1)*nObs,1+seq(nModpar)] <- 
                        fdHess(parModel, .model.solve.fdHess.interface,
                               model=model, t=time[it], iObs=iObs, 
                               iArms)$gradient
                }
            }
        }
    }
    res
}

model.sensJac <- function(model, time=model$time) {
    .model.solve.jacUser.interface <- function(tim, y, model) {
        tmp <- model$PFIMmodel(tim, y[seq(model$nIC)], model$parModel, model$parArmsNow)
                                        # Normal right hand side
        dy <- tmp[[1]]
        ## Partial derivatives...
        ## dxi/dpj -> d/dt(dxi/dpj) = d/dpj(dxi/dt) =
        ##                d/dpj(fi) = Sum[@fi/@xk*dxk/dpj] + @fi/@pj
        JacState <- tmp[[3]]; JacParam <- tmp[[4]] 
        pDer <- y[model$nIC+seq(model$nIC*model$nModpar)]; dim(pDer) <- c(model$nIC,model$nModpar)
        dpDer <- JacState %*% pDer + JacParam
        
                                        # If observations are not just the state vectors
        Obs <- NULL
        dObs <- NULL
        if (length(tmp)>4) {
            Obs <- tmp[[2]]      # Solution vector for observations
            dObs <- tmp[[6]] %*% pDer + tmp[[5]] # Partial derivatives of observations
        }
        
        list(c(dy, dpDer), Obs, dObs)
    }
    
                                        # Define user-friendly parameters
    nModpar <- model$nModpar
    nObs <- model$nObs
    nArms <- model$nArms
    parModel <- model$parModel
    parArms <- model$parArms

    rtol <- atol <- sqrt(.Machine$double.eps) # default: 1e-6
    hmax <- model$Hmax
                                        # Set initial time to 0 if min(time) > 0
    setInitialTime <- FALSE
    if (time[1] > 0) {
        time <- c(0, time)
        setInitialTime <- TRUE
    }
	
    ## Define general model parameters (for use in initial conditions)
    for (iModpar in seq(nModpar)) assign(model$parModelName[iModpar], parModel[iModpar])
	
    res <- array(dim=c(length(time), nObs*nArms, 1+nModpar))
    for (iArms in seq(nArms)) {
        ## Define arm/group parameters (for use in initial conditions)
        for (iArmpar in seq(model$nArmspar))
            assign(model$parArmsName[iArmpar], parArms[iArmpar,iArms])
                                        # Derive initial conditions
        IC <- eval(model$parArmsIC[[iArms]])
        model$nIC <- length(IC)
                                        # Include ICs for partial derivatives
        if (length(JacIC<-eval(model$JacIC[[iArms]])) == 0)
            stop(paste("Analytical sensitivities with ODES selected: ",
                       "Supposed to provide partial differentials, ",
                       "but missing model$JacIC, derivatives of initial conditions", collapse='\n'))
        IC <- c(IC, JacIC)
        
                                        # If dosing event entry exists, then add to lsoda
        if (!is.null(model$dosing)) {
            names(IC)[1] <- "v1"
            eventdat <- model$dosing[[iArms]]
        } else eventdat <- NULL
        
                                        # Solve ODE system (with analytical partial derivatives)
        model$parArmsNow <- model$parArms[[iArms]]
        tmp <- lsoda(IC, time, .model.solve.jacUser.interface, model, rtol=rtol, atol=atol, hmax=hmax, events = list(data = eventdat))[,-1]
        
        ## If observations are not the state variables but derived variables
        if (dim(tmp)[2] != nObs*(1+nModpar)) {
            tmp <- tmp[, model$nIC*(1+nModpar)+seq(nObs*(1+nModpar))]
        }
        
        dim(tmp) <- c(length(time), nObs, 1+nModpar)
        res[, seq(nObs) + (iArms-1)*nObs,] <- tmp
    }
    ## If initial time was set to 0 (without having a measurement
    ## point there, we have to remove it again)
    if (setInitialTime) res <- res[-1,,,drop=F]
    
    res
}

