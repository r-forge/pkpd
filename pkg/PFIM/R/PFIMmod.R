##' Expand information in the PFIMmod structure
##'
##' 
##' @title Create a PFIM model structure
##' @param PFIMmodel either a list containing the elements of the model or a function or an expression
##' @param Type character - one of 'AE' (algebraic equations) or 'DE' differential equations
##' @param parModelName character vector of model parameter names 
##' @param parModel numeric - model parameter values
##' @param parModelVar numeric - model parameter variances
##' @param parModelVarType character - one of 'add' or 'exp'
##' @param parObsName - character vector of observation/response names
##' @param parObsErr error amount and structure related to observations
##' @param parObsTimes list of possible sampling times per observation
##' @param parArms arm/group-specific parameters
##' @param parArmsName names of arm/group-specific parameters
##' @param ArmsName labels of the arms
##' @param TimeName label for the time axis
##' @param tRange numeric vector of length 2 giving the range of possible times
##' @param parArmsIC initial conditions (for ODEs only)
##' @return list of model specification
PFIMmod <- function(PFIMmodel, Type, parModelName, parModel, parModelVar,
                    parModelVarType, parObsName, parObsErr, parObsTimes,
                    parArms, parArmsName, ArmsName, TimeName, tRange, parArmsIC) {
    if (is.list(PFIMmodel)) {
        model <- PFIMmodel
    } else {
                                        # Collate given input parameters
        model <-list(PFIMmodel=PFIMmodel, Type=Type, parModelName=parModelName,
                     parModel=parModel, parModelVar=parModelVar,
                     parModelVarType=parModelVarType, parObsName=parObsName,
                     parObsErr=parObsErr, parObsTimes=parObsTimes,
                     parArms=parArms, parArmsName=parArmsName,
                     ArmsName=ArmsName, TimeName=TimeName, tRange=tRange,
                     parArmsIC=if(Type == "DE") parArmsIC else NULL)
    }
    stopifnot(model$Type %in% c("AE", "DE"))
                                        # Number of parameters and observations
    model$nModpar <- length(model$parModelName)
    model$nObs <- length(model$parObsName)
    model$nArmspar <- length(model$parArmsName)
    model$nArms <- length(model$ArmsName)
                                        # Add info for graphical logarithmic output
                                        # Assume that first obs most often PK
    if (is.null(model$parObsLogY))
        model$parObsLogY <- c(TRUE, rep(FALSE, model$nObs - 1L))
    if (is.null(model$parObsLogYMin))
        model$parObsLogYMin <- rep(1e-4, model$nObs)
                                        # Time components (time, observations, arms)
    model$time <- sort(unique(c(unlist(model$parObsTimes), # automatically handles null model$dosing
                                vapply(model$dosing, function(mm)mm[,1], 0, USE.NAMES=FALSE))))
    model$nTime <- length(model$time)   # probably not necessary
    model$fullTime <- sort(unique(c(seq(min(c(model$tRange[1], model$time[1])),
                                        max(model$tRange[2], model$time[2]),
                                        length.out=100), model$time)))

    tmp <- data.frame(time=rep(model$time, times=model$nObs),
                      iObs=rep(1:model$nObs, each=model$nTime))
    if (!is.null(model$parObsTimes$PFIM)) { # PFIM - sorted by observations
        for (iObs in 1:(length(model$parObsTimes)-1)) {
            for (iArm in 1:length(model$parObsTimes[[iObs]])) {
                tmp[tmp$iObs==iObs, 2+iArm] <-
                    (model$time %in% model$parObsTimes[[iObs]][[iArm]])
            }
        }
        if (model$nArms>1 & ncol(tmp)==3) {
            ## If for all observations only time points for first arm
            ## are given - extend/copy to all arms
            tmp[,2+2:model$nArms] <- tmp[,3]
        }
    } else { # my structure - sorted by arm
        stopifnot(length(model$parObsTimes) >= model$nArms)
        for (iArm in 1:length(model$parObsTimes)) {
            if (is.list(model$parObsTimes[[iArm]])) {
                for (iObs in 1:length(model$parObsTimes[[iArm]])) {
                    tmp[tmp$iObs==iObs, 2+iArm] <-
                        (model$time %in% model$parObsTimes[[iArm]][[iObs]])
                }
            } else {
                ## If for an arms only one set of time points are
                ## given - use/copy for all observations 
                tmp[, 2+iArm] <- (model$time %in% model$parObsTimes[[iArm]])
            }
        }
    }
    model$timeMap <- tmp

    ## Create matrices out of list objects (for being user friendly and avoid confusion!)
    if (is.list(model$parObsErr)) {
        stopifnot((nObsErr<-length(unlist(model$parObsErr))) %% 2 != 1)
        model$parObsErr <- matrix(unlist(model$parObsErr), nrow=2, ncol=nObsErr/2)
    }
    if (is.list(model$parArms)) {
        nParArms <- length(model$parArmsName)
        nArms <- length(model$ArmsName)
        stopifnot(length(unlist(model$parArms)) == nParArms*nArms)
        model$parArms <- matrix(unlist(model$parArms), nrow=nParArms, ncol=nArms, byrow=TRUE)
    }
                                        # For DE ensure that respective libraries are loaded
    if (model$Type == 'DE')	{
        require(deSolve)
        if (!is.null(model$dosing))
            if (!is.data.frame(model$dosing)) {
                ## If dosing is data.frame we assume that it is done according to lsoda needs...
                ## Otherwise we implement all doses into the first compartment
                eventdat <- list()
                for (iArms in seq(nArms)) {
                    eventdat[[iArms]] <- data.frame(var = "v1",
                                                    time = model$dosing[[iArms]][,1],
                                                    value = model$dosing[[iArms]][,2],
                                                    method = "add")
                }
                model$dosing <- eventdat
            }
    }
                                        # PFIM backward compatibility
    if (model$parModelVarType %in% c(1,2)) {
        model$parModelVarType <- c('add','exp')[model$parModelVarType]
        warning(paste("Random effect model should be either 'add' or 'exp'.",
                      "('1' and '2' currently allowed for backward compatibility)",
                      sep="\n"))
    }

    class(model) <- c(paste("PFIMmod", model$Type, sep=""), "PFIMmod")
    check(model)
}

check <- function(object) UseMethod("check")

##' Checks model input and output consistency
##' Runs the model file once for getting number of derivatives, etc.
##'
##' Checks run include:
##' Run model once
##' Check input consistency
##' Check output consistency with specifications
##' @title Check an object for consistency
##' @param object 
##' @return object
check.PFIMmod <- function(object) {
    with(object, stopifnot(             # Model parameter specifications
                           length(parModelName) == length(parModel),
                           length(parModelName) == length(parModelVar),
                           all(parModel != 0),
                           parModelVarType %in% c('exp','add'),
                           all(parModelVar >= 0),
                                        # Observation specifications
                           length(parObsName) ==ncol(parObsErr),
                           2L == nrow(parObsErr),
                           length(parObsLogY) == nObs,
                           length(parObsLogYMin) == nObs,
                                        # Arms/groups specifications
                           length(parArmsName) == nrow(parArms),
                           length(ArmsName) == ncol(parArms),
                                        # General model type, structure, etc.
                           Type %in% c('DE','AE'))
         )
    object
}

check.PFIMmodAE <- function(object) {
                                        # Check length of observations with number of equations
    with(object, if (is.expression(PFIMmodel)) {
        stopifnot(nObs == length(eval(PFIMmodel)))
        NextMethod()
    })
    stopifnot(is.function(object$PFIMmodel))
    with(object,
         if ("model" %in% names(formals(PFIMmodel)))
     {
         stopifnot(nObs ==
                   length(PFIMmodel(0,parModel,parArms[,1], object)))
     } else {
         parMod <- c(parModel)
         names(parMod) <- parModelName
         lst <- as.list(parMod)
         lst[[parArmsName]] <- as.vector(parArms)
         lst[[TimeName]] <- parObsTimes[[1]][[1]]
         stopifnot(length(lst[[TimeName]]) == length(do.call(PFIMmodel, lst)))
     })
    NextMethod()
}
check.PFIMmodDE <- function(object) {
                                        # Do something
    if (!is.null(object$dosing)) {
        if (object$nArms != length(object$dosing))
            stop("Number of entries in object$dosing must coincide with number of arms specified.")
        if (!is.data.frame(object$dosing[[1]]))
            for (iArms in seq(object$nArms)) {
                if (ncol(object$dosing[[iArms]]) != 2)
                    stop("Columns of each entry in object$dosing must be for time and respective dose, i.e., 2 columns.")
            }
    }
                                        # Check validity of ICs
    ## for (i in 1:p) {assign(parameters[i],theta[i]) }
    ## cond<-eval(condinit[kit])	
    NextMethod()
}

solve.PFIMmodAE <- function(a, b, ...) { # this is probably not a good generic to co-opt
                                        # remap the names from the generic
    model <- a
    tim <- if(missing(b)) model$time else b
    calcSens <- list(...)$calcSens
    if (is.null(calcSens)) calcSens <- FALSE

    nModpar <- model$nModpar
    parModel <- model$parModel
    nObs <- model$nObs
    parArms <- model$parArms
    nArms <- model$nArms
    nArmspar <- model$nArmspar

    if (calcSens) {
        res <- array(dim=c(length(tim), nObs*nArms, 1+nModpar))
    } else {
        res <- array(dim=c(length(tim), nObs*nArms, 1))
    }
### FIXME: This is better done with an environment
                                        # Assign model parameter values
    for (iModpar in seq(nModpar)) assign(model$parModelName[iModpar], parModel[iModpar])
                                        # For PFIM-like expression as model
    if (is.expression(model$PFIMmodel)) {
        t <- tim
        modelEQ <- model$PFIMmodel
        for (iArms in seq(nArms)) {
            for (iArmpar in seq(nArmspar))
                assign(model$parArmsName[iArmpar], parArms[iArmpar,iArms])
            for (iObs in seq(nObs)) {
                res[, iObs + (iArms-1)*nObs, 1] <- eval(modelEQ[iObs])
            }
        }

        if (calcSens==TRUE) {
            for (iArms in seq(nArms)) {
                for (iArmpar in seq(nArmspar))
                    assign(model$parArmsName[iArmpar], parArms[iArmpar,iArms])
                for (iModpar in seq(nModpar)) {
                    for (iObs in seq(nObs)) {
                        modelEQ <- D(model$PFIMmodel[iObs], model$parModelName[iModpar])
                        res[, iObs + (iArms-1)*nObs, 1+iModpar] <- eval(modelEQ)
                    }
                }
            }
        }
        return(res)
    }
                                        # For model <- function(tim, parModel, parArms, opt)
    if (is.function(model$PFIMmodel)) {
        ## If sensitivities should be calculated then do this by numerical derivatives
        if (calcSens==TRUE)
            return(model.sensNum(model, tim))
        ## This provides the solution of the system by function call
        for (iArms in seq(nArms)) {
            tmp <- try(withCallingHandlers(model$PFIMmodel(tim, parModel, parArms[,iArms], model),
                                           warning = stop()), silent=TRUE)
            if (is.numeric(tmp) & length(tmp)==length(tim)*nObs) {
                ## If the function is 'vectorized' it immediately
                ## prints out all observations at all times
                res[,(iArms-1)*nObs + seq(nObs),1] <- tmp
            } else {
                ## If the function is not as flexible we have to cycle through all time points
                res[,(iArms-1)*nObs + seq(nObs),1] <-
                    t(sapply(1:length(tim),
                             function(tidx, parModel, parA, model)
                             model$PFIMmodel(tim[tidx], parModel, parA, model),
                             parModel=parModel,
                             parA=parArms[,iArms],
                             model=model))
            }
        }
        return(res)
    }
### FIXME: Why not call stop instead of warning?  I would say this is
### an error.
    warning("Could not determine model type (expression, function, ODE,...) properly")
    NULL
}

solve.PFIMmodDE <- function(a, b, ...) {
                                        # remap the names from the generic
    model <- a
    tim <- if(missing(b)) model$time else b
    calcSens <- list(...)$calcSens
    if (is.null(calcSens)) calcSens <- FALSE

    modelDE <- model$PFIMmodel
    nModpar <- model$nModpar
    nObs <- model$nObs
    parArms <- model$parArms
    nArms <- model$nArms
    nArmspar <- model$nArmspar
    parModel <- model$parModel

    atol <- rtol <- sqrt(.Machine$double.eps) # default: 1e-6
    hmax <- model$Hmax
                                        # Set initial time to 0 if min(time) > 0
    setInitialTime <- FALSE
    if (time[1] > 0) {
        time <- c(0, time)
        setInitialTime <- TRUE
    }

    ## Define general model parameters (for possible use in initial conditions)
    for (iModpar in seq(nModpar)) assign(model$parModelName[iModpar], parModel[iModpar])

    res <- array(dim=c(length(time), nObs*nArms, 1))
    for (iArms in seq(nArms)) {
        ## Define arm/group parameters (for possible use in initial conditions)
        for (iArmpar in seq(nArmspar)) assign(model$parArmsName[iArmpar], parArms[iArmpar,iArms])
        IC <- eval(model$parArmsIC[[iArms]])
        obsIdx <- 1 + length(IC) + seq(nObs)
                                        # If dosing event entry exists, then add to lsoda
        if (!is.null(model$dosing)) {
            names(IC)[1] <- "v1"
            eventdat <- model$dosing[[iArms]]
        } else eventdat <- NULL
        res[, seq(nObs) + (iArms-1)*nObs, 1] <-
            lsoda(IC, time, modelDE, parModel, rtol=rtol, atol=atol, hmax=hmax,
                  events=list(data=eventdat))[, obsIdx]
    }

    ## If initial time was set to 0 (without having a measurement
    ## point there, we have to remove it again)
    if (setInitialTime) {
        res <- res[-1, , , drop=F]
    }
    res
}
