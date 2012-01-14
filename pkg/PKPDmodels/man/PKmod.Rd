\name{PKmod}
\alias{PKmod}
\title{PK models with linear elimination}
\usage{
  PKmod(admin = c("bolus", "infusion", "oral"),
    dosage = c("sd", "md", "ss"),
    cpt = c("one", "two", "three"), subst = list(),
    hessian = FALSE)
}
\arguments{
  \item{admin}{form of administration of the drug, one of
  \code{"bolus"}, \code{"infusion"} or \code{"oral"}.
  Defaults to}

  \item{dosage}{form of dosage, one of \code{"sd"} (single
  dose), \code{"md"} (multiple, equally-spaced doses) and
  \code{"ss"} (steady-state).  Defaults to \code{"sd"}.}

  \item{cpt}{number of model compartments, one of
  \code{"one"}, \code{"two"} or \code{"three"}.  Defaults
  to \code{"one"}.}

  \item{subst}{a list of formulas of substitutions to
  perform}

  \item{hessian}{a logical value indicating whether the
  second derivatives should be calculated and incorporated
  in the return value.}
}
\value{
  a byte-compiled model function with gradient evaluation
}
\description{
  Create a model function with gradient evaluation (and,
  optionally, Hessian evaluation) for a model according to
  the number of compartments, the form of administration
  and dosage of the drug after performing any substitutions
  given.
}
\details{
  The substitutions are given as a list of formulas, such
  as \code{list(k ~ Cl/V, Cl ~ exp(lCl), V ~ exp(lV))}.
  They are applied left to right.
}
\examples{
## return a function with substitutions
PKmod("bolus", "sd", "one", list(k ~ Cl/V, Cl ~ exp(lCl), V ~ exp(lV)))
}

