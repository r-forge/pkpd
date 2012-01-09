\name{sensitivity}
\alias{sensitivity}
\title{Evaluate sensitivities (gradients)}
\usage{sensitivity(model, time=model$time, useHess=FALSE, ...)}
\description{Evaluate the sensitivities (gradients) of a PFIM model}
\details{Evaluate the sensitivities of a PFIM model object with respect to
the model parameters.}
\value{a three-dimensional array of model function values and sensitivities}
\arguments{\item{model}{a PFIM model object}
\item{time}{an optional numeric vector of times at which to evaluate}
\item{useHess}{logical - whether to use the \code{fdHess} function from
the \pkg{nlme} package or the \code{gradient} function from the
\pkg{numDeriv} package.}
\item{\dots}{optional arguments, none used at present}}
