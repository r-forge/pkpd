\name{npde-package}
\alias{npde-package}
\docType{package}
\title{Normalised prediction distribution errors for nonlinear mixed-effect models}
\description{Routines to compute normalised prediction distribution errors, a metric designed to evaluate non-linear mixed effect models such as those used in pharmacokinetics and pharmacodynamics}
\details{
\tabular{ll}{
Package: \tab npde\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-02-14\cr
License: \tab GPL version 2 or later\cr
}
  See the documentation for npde or autonpde for details
}
\author{
Emmanuelle Comets, Karl Brendel and France Mentre

Maintainer: Emmanuelle Comets <emmanuelle.comets@bichat.inserm.fr>
}
\references{K. Brendel, E. Comets, C. Laffont, C. Laveille, and F. Mentre. Metrics for external model evaluation with an application to the population pharmacokinetics of gliclazide. \emph{Pharmaceutical Research}, 23:2036--49, 2006.}
\keyword{models}
\examples{
data(theopp)
data(simtheopp)

# Calling autonpde with dataframes

# x<-autonpde(theopp,simtheopp,1,3,4,boolsave=FALSE)
# x$npde

# Calling autonpde with names of files to be read from disk

# write.table(theopp,"theopp.tab",quote=FALSE,row.names=FALSE)
# write.table(simtheopp,"simtheopp.tab",quote=FALSE,row.names=FALSE)
# x<-autonpde(namobs="theopp.tab", namsim="simtheopp.tab", iid = 1, ix = 3, iy = 4, boolsave = FALSE)
# x$npde

}
