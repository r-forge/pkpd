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

##' Return a formula for the one-compartment PK model with linear
##' elimination according to the administration form and the dosage pattern.
##' 
##' @title Expressions for one-compartment PK models with linear elimination
##'
##' @param admin form of administration of the drug, one of
##'    \code{"bolus"}, \code{"infusion"} or \code{"oral"}.  Defaults to
##     \code{"bolus"}. 
##' @param dosage form of dosage, one of \code{"sd"} (single dose),
##'    \code{"md"} (multiple, equally-spaced doses) and \code{"ss"}
##'    (steady-state).  Defaults to \code{"sd"}.
##' @param subst a list of formulas of substitutions to perform
##' @return a formula
##' @examples
##' ## single-dose oral administration
##' PK1expr("oral", "sd")
PK1expr <- function(admin=c("bolus", "infusion", "oral"),
                    dosage=c("sd", "md", "ss"), subst=list()) {
    frm <- list(bolus =
                list(sd = ~dose * exp(-k * t) / V,
                     md = ~(dose/V) * ((1-exp(-N*k*tau))/(1-exp(-k*tau))) * exp(-k*(t-(N-1)*tau)),
                     ss = ~(dose/V)/(1-exp(-k*tau))*(exp(-k*(t-(TimeSS))))),
                infusion =
                list(),
                oral =
                list(sd = ~(dose/V) * (ka/(ka-k)) * (exp(-k*t)-exp(-ka*t)),
                     md = ~(dose/V) * (ka/(ka-k)) * (exp(-k*(t-(N-1)*tau))*(1-exp(-N*k*tau)) /
                                                     (1-exp(-k*tau))-exp(-ka*(t-(N-1)*tau)) *
                                                     (1-exp(-N*ka*tau))/(1-exp(-ka*tau))),
                     ss = ~(dose/V) * (ka/(ka-k)) * (exp(-k*(t-TimeSS))/(1-exp(-k*tau))-
                                                     exp(-ka*(t-TimeSS))/(1-exp(-ka*tau))))
                )[[match.arg(admin)]][[match.arg(dosage)]]
    for (i in seq_along(subst)) {
        stopifnot(class(subfrm <- eval(subst[[i]])) == "formula",
                  is.name(subfrm[[2]]))
        frm <- subexpr(frm, subfrm[[2]], as.expression(subfrm[[3]]))
    }
    frm
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
##' @param admin form of administration of the drug, one of
##'    \code{"bolus"}, \code{"infusion"} or \code{"oral"}.  Defaults to
##     \code{"bolus"}. 
##' @param dosage form of dosage, one of \code{"sd"} (single dose),
##'    \code{"md"} (multiple, equally-spaced doses) and \code{"ss"}
##'    (steady-state).  Defaults to \code{"sd"}.
##' @param subst a list of formulas of substitutions to perform
##' @return a byte-compiled model function with gradient evaluation 
##'
##' @examples
##' ## return a function with substitutions
##' PK1cmpt("bolus", "sd", list(k ~ Cl/V, Cl ~ exp(lCl), V ~ exp(lV)))
##'
PK1cmpt <- function(admin=c("bolus", "infusion", "oral"),
                    dosage=c("sd", "md", "ss"),
                    subst=list()) {
    frm <- PK1expr(admin, dosage, subst)
    covariates <- c("dose", "t",
                    list(sd=character(0), md=c("N", "tau"), ss=c("TimeSS", "tau"))[[dosage]])
    pnms <- setdiff(all.vars(frm), covariates)
    cmpfun(deriv(frm, pnms, c(covariates, pnms)))
}
