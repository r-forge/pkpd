\name{PK1cmpt}
\alias{PK1cmpt}
\title{One-compartment PK models with linear elimination}
\usage{PK1cmpt(admin, subst=list(), formulaOnly=FALSE)}
\description{Create a model function with gradient evaluation for a
one-compartment model according to the form of administration of
the drug after performing any substitutions given.}
\details{The substitutions are given as a list of formulas, such as 
list(k ~ Cl/V, Cl ~ exp(lCl), V ~ exp(lV)).  They are applied left
to right.}
\value{a model function with gradient evaluation unless
\code{formulaOnly} is \code{TRUE}}
\arguments{\item{admin}{form of administration of the drug, a character value
used to index the list \code{forms} defined in the function.}
\item{subst}{a list of formulas of substitutions to perform}
\item{formulaOnly}{logical - when TRUE the formula only is returned}}
\examples{PK1cmpt("bolus", formulaOnly=TRUE)  # simplest usage
## return a function with substitutions
PK1cmpt("bolus", list(k ~ Cl/V, Cl ~ exp(lCl), V ~ exp(lV)))}
