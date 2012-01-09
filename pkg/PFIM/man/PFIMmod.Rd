\name{PFIMmod}
\alias{PFIMmod}
\title{Create a PFIM model structure}
\usage{PFIMmod(PFIMmodel, Type, parModelName, parModel, parModelVar,
    parModelVarType, parObsName, parObsErr, parObsTimes, parArms,
    parArmsName, ArmsName, TimeName, tRange, parArmsIC)}
\description{Expand information in the PFIMmod structure}
\value{list of model specification}
\arguments{\item{PFIMmodel}{either a list containing the elements of the model or a function or an expression}
\item{Type}{character - one of 'AE' (algebraic equations) or 'DE' differential equations}
\item{parModelName}{character vector of model parameter names}
\item{parModel}{numeric - model parameter values}
\item{parModelVar}{numeric - model parameter variances}
\item{parModelVarType}{character - one of 'add' or 'exp'}
\item{parObsName}{- character vector of observation/response names}
\item{parObsErr}{error amount and structure related to observations}
\item{parObsTimes}{list of possible sampling times per observation}
\item{parArms}{arm/group-specific parameters}
\item{parArmsName}{names of arm/group-specific parameters}
\item{ArmsName}{labels of the arms}
\item{TimeName}{label for the time axis}
\item{tRange}{numeric vector of length 2 giving the range of possible times}
\item{parArmsIC}{initial conditions (for ODEs only)}}
