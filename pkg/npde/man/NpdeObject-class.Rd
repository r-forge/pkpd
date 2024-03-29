\docType{class}
\name{NpdeObject-class}
\alias{npde.graphs,NpdeObject}
\alias{npde.main,NpdeObject}
\alias{NpdeObject}
\alias{NpdeObject-class}
\alias{NpdeObject-class,}
\alias{[<-,NpdeObject-method}
\alias{[,NpdeObject-method}
\alias{npde.save,NpdeObject}
\alias{plot,NpdeObject}
\alias{print,NpdeObject-method}
\alias{showall,NpdeObject-method}
\alias{show,NpdeObject-method}
\alias{summary,NpdeObject-method}
\alias{test,NpdeObject-method}
\title{Class "NpdeObject"}
\description{
  An object of class NpdeObject
}
\section{Objects from the Class}{
  NpdeObject objects are typically created by calls to
  \code{\link{npde}} or \code{\link{autonpde}}. They
  contain the following slots:

  \describe{ \item{data}{an object of class NpdeData,
  containing the observed data} \item{sim.data}{an object
  of class NpdeSimData, containing the simulated data}
  \item{res}{an object of class NpdeRes, containing the
  results} \item{options}{a list of options} \item{prefs}{a
  list of graphical preferences for the plots} }
}

\section{Methods}{
  \describe{ \item{print(x):}{Prints a summary of object}
  \item{show(x):}{Prints a short summary of object}
  \item{showall(x):}{Prints a detailed summary of object}
  \item{plot(x):}{Diagnostic and other plots. More details
  can be found in \code{\link{plot.NpdeObject}}}
  \item{summary(x):}{Returns a summary of object x in list
  format} \item{gof.test(x, which="npde", parametric=TRUE,
  ...):}{Returns goodness-of-fit tests}
  \item{set.npde.plotoptions(x):}{Sets options for graphs
  (internal method used in plots)} }
}
\examples{
methods(class="NpdeObject")

showClass("NpdeObject")
}
\seealso{
  \code{\link{npde}}, \code{\link{autonpde}},
  \code{\link{NpdeData}}, \code{\link{NpdeSimData}},
  \code{\link{NpdeRes}}, \code{\link{gof.test}}
}
\keyword{classes}

