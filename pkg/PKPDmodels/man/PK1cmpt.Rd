\name{PK1cmpt}
\alias{PK1cmpt}
\title{One-compartment PK models with linear elimination}
\usage{PK1cmpt(admin=c("bolus", "infusion", "oral"),
    dosage=c("sd", "md", "ss"), subst=list())}
\description{Create a model function with gradient evaluation for a
one-compartment model according to the form of administration of
the drug after performing any substitutions given.}
\details{The substitutions are given as a list of formulas, such as 
\code{list(k ~ Cl/V, Cl ~ exp(lCl), V ~ exp(lV))}.  They are applied left
to right.}
\value{a byte-compiled model function with gradient evaluation}
\arguments{\item{admin}{form of administration of the drug, one of
\code{"bolus"}, \code{"infusion"} or \code{"oral"}.  Defaults to}
\item{dosage}{form of dosage, one of \code{"sd"} (single dose),
\code{"md"} (multiple, equally-spaced doses) and \code{"ss"}
(steady-state).  Defaults to \code{"sd"}.}
\item{subst}{a list of formulas of substitutions to perform}}
\examples{## return a function with substitutions
PK1cmpt("bolus", "sd", list(k ~ Cl/V, Cl ~ exp(lCl), V ~ exp(lV)))}
