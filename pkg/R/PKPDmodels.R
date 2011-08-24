##' Substitute the expression \code{sub} for the name \code{nm} in
##' \code{expr} by walking the tree.
##'
##' @title "Sub[stitute] expr[ession]"
##'
##' @param expr an expression
##' @param nm a name for which to substitute
##' @param sub the expression to substitute for name nm
##'
##' @return the expression with all occurrences of nm replaced by sub
##' @note this function is called recursively
subexpr <- function(expr, nm, sub) 
{
    if (length(expr) == 1) {
        if (is.name(expr) && expr == nm) return(sub[[1]])
        return(expr)
    }
    for (j in 2:length(expr)) expr[[j]] <- subexpr(expr[[j]], nm, sub)
    expr
}
    
##' Create a model function with gradient evaluation for a
##' one-compartment model according to the form of administration of
##' the drug after performing any substitutions given. 
##' 
##'
##' The substitutions are given as a list of formulas, such as 
##' \code{list(k ~ Cl/V, Cl ~ exp(lCl), V ~ exp(lV))}.  They are applied left
##' to right.
##' 
##' @title One-compartment PK models with linear elimination
##'
##' @param admin form of administration of the drug, a character value
##'    used to index the list \code{forms} defined in the function.
##' @param subst a list of formulas of substitutions to perform
##' @param formulaOnly logical - when TRUE the formula only is returned
##'
##' @return a model function with gradient evaluation unless
##'    \code{formulaOnly} is \code{TRUE}
##'
##' @examples
##' PK1cmpt("bolus", formulaOnly=TRUE)  # simplest usage
##'    ## return a function with substitutions
##' PK1cmpt("bolus", list(k ~ Cl/V, Cl ~ exp(lCl), V ~ exp(lV)))
##'
PK1cmpt <- function(admin, subst=list(), formulaOnly=FALSE) {
    stopifnot(is.character(admin), is.list(subst))
    forms <- list(bolus = ~dose * exp(-k * t) / V,
                  mbolus = ~(dose/V) * ((1-exp(-N*k*tau))/(1-exp(-k*tau))) * exp(-k*(t-(N-1)*tau)),
                  ssbolus = ~(dose/V)/(1-exp(-k*tau))*(exp(-k*(t-(TimeSS)))),
                  oral = ~(dose/V) * (ka/(ka-k)) * (exp(-k*t)-exp(-ka*t)), # single oral dose
                  moral = ~ (dose/V) * (ka/(ka-k)) * # multiple oral doses
                  (exp(-k*(t-(N-1)*tau))*(1-exp(-N*k*tau)) /
                   (1-exp(-k*tau))-exp(-ka*(t-(N-1)*tau))*(1-exp(-N*ka*tau))/(1-exp(-ka*tau))),
                  ssoral = ~ (dose/V) * (ka/(ka-k)) * (exp(-k*(t-TimeSS))/(1-exp(-k*tau))-
                                                       exp(-ka*(t-TimeSS))/(1-exp(-ka*tau))))
    frm <- forms[[admin]]
    for (i in seq_along(subst)) {
        stopifnot(class(subfrm <- eval(subst[[i]])) == "formula",
                  is.name(subfrm[[2]]))
        frm <- subexpr(frm, subfrm[[2]], as.expression(subfrm[[3]]))
    }
    if (formulaOnly) return(frm)
    pnms <- setdiff(all.vars(frm), c("dose", "t", "N", "tau", "TimeSS"))
    deriv(frm, pnms, c("dose", "t", pnms))
}

