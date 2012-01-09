##' Read a PFIM3.2.2 model specification and create the internal representation
##'
##' 
##' @title Read a PFIM3.2.2 model specification
##' @param dir directory containing files named 'stdin.r' and 'model.r'
##' @param model_created logical is there another file called 'model_created.r'
##' @return a PFIMmod object
readPFIM <- function(dir, model_created=FALSE) {
                                        # Read in default file
    source(file.path(dir, "stdin.r"), local=TRUE)
    source(file.path(dir, file.model), local=TRUE)
                                        # Current restrictions!
    if (covariate.model || covariate_occ.model || n_occ > 1)
        stop('Covariate model and IOV not implemented, yet.')
    stopifnot(run == 'EVAL', option == 1, nr == length(names.datay))
    Obsnms <- LETTERS[seq(as.integer(nr))]
    ObsErr <- lapply(Obsnms, function(nm)
                     c(get(paste("sig.inter", nm, sep='')),
                       get(paste("sig.slope", nm, sep=''))))
    ObsTimes <- lapply(Obsnms, function(nm) get(paste("prot", nm, sep='')))
    Type <- c(AF="AE", DE="DE")[modelform]

    if (exists('dose') & Type == 'AE') {
        parArms <- as.list(dose)
        parArmsName <- 'dose'
        ArmsName <- as.list(paste('Dose',dose))
    } else if (Type == 'DE') {
        nArms <- length(condinit)
        parArms <- rep(list(c(0)), nArms)
        parArmsName <- c('dummy')
        ArmsName <- as.list(paste('Output',seq(nArms)))
    } else {
        parArms <- list(c(0))
        parArmsName <- c('dummy')
        ArmsName <- list('Output')
    }
                                        # construct initial value of the list
    PFIMmod(switch(Type, AE=form, DE=formED), Type, parameters, beta, diag(omega),
            Trand, names.datay, ObsErr, ObsTimes, parArms, parArmsName, ArmsName,
            names.datax[1], c(0, max(unlist(ObsTimes))), NULL)
}
