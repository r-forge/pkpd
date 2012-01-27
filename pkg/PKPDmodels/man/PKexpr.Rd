\name{PKexpr}
\alias{PKexpr}
\title{Expressions for PK models with linear elimination}
\usage{PKexpr(admin=c("bolus", "infusion", "oral"), dosage=c("sd", "md", "ss"),
    cpt=c("one", "two", "three"), subst=list())}
\description{Return a formula for the PK model with linear
elimination according to the number of compartments, the administration form and the dosage pattern.}
\value{a formula}
\arguments{\item{admin}{form of administration of the drug, one of
\code{"bolus"}, \code{"infusion"} or \code{"oral"}.  Defaults to}
\item{dosage}{form of dosage, one of \code{"sd"} (single dose),
\code{"md"} (multiple, equally-spaced doses) and \code{"ss"}
(steady-state).  Defaults to \code{"sd"}.}
\item{cpt}{number of model compartments, one of \code{"one"},
\code{"two"} or \code{"three"}.  Defaults to \code{"one"}.}
\item{subst}{a list of formulas of substitutions to perform}}
\examples{## single-dose oral administration
PKexpr("oral", "sd", "one")}
