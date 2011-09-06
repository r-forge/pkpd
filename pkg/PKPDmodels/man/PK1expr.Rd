\name{PK1expr}
\alias{PK1expr}
\title{Expressions for one-compartment PK models with linear elimination}
\usage{PK1expr(admin=c("bolus", "infusion", "oral"), dosage=c("sd", "md",
    "ss"), subst=list())}
\description{Return a formula for the one-compartment PK model with linear
elimination according to the administration form and the dosage pattern.}
\value{a formula}
\arguments{\item{admin}{form of administration of the drug, one of
\code{"bolus"}, \code{"infusion"} or \code{"oral"}.  Defaults to}
\item{dosage}{form of dosage, one of \code{"sd"} (single dose),
\code{"md"} (multiple, equally-spaced doses) and \code{"ss"}
(steady-state).  Defaults to \code{"sd"}.}
\item{subst}{a list of formulas of substitutions to perform}}
\examples{## single-dose oral administration
PK1expr("oral", "sd")}
