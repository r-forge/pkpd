\name{PKexpr}
\alias{PKexpr}
\title{Expressions for PK models with linear elimination}
\usage{
  PKexpr(admin = c("bolus", "infusion", "oral"),
    dosage = c("sd", "md", "ss"), subst = list(), cpt = 1L)
}
\arguments{
  \item{admin}{form of administration of the drug, one of
  \code{"bolus"}, \code{"infusion"} or \code{"oral"}.
  Defaults to \code{"bolus"}.}

  \item{dosage}{form of dosage, one of \code{"sd"} (single
  dose), \code{"md"} (multiple, equally-spaced doses) and
  \code{"ss"} (steady-state).  Defaults to \code{"sd"}.}

  \item{subst}{a list of formulas of substitutions to
  perform}

  \item{cpt}{scalar integer - the number of model
  compartments.}
}
\value{
  a formula
}
\description{
  Return a formula for the PK model with linear elimination
  according to the number of compartments, the
  administration form and the dosage pattern.
}
\examples{
## single-dose oral administration
PKexpr("oral", "sd")
}

