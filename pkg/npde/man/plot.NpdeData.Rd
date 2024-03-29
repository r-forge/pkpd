\name{plot.NpdeData}
\alias{plot.NpdeData}
\title{Plots a NpdeData object}
\usage{
  \method{plot}{NpdeData} (x, y, ...)
}
\arguments{
  \item{x}{a NpdeData object}

  \item{y}{unused, here for compatibility with the base
  plot function}

  \item{\dots}{additional graphical parameters to be passed
  on to the plot}
}
\description{
  Plots the data in a NpdeData object
}
\details{
  The default plot is a spaghetti plot of all the data,
  with a line joining the observations for each subject. If
  censored data is present, it is shown with a different
  symbol and colour.
}
\examples{
data(theopp)

x<-npdeData(theopp,name.group="ID",name.predictor="Time",name.response="Conc",
name.covariates=c("Wt"),units=list(x="hr",y="mg/L",covariates="kg"))
plot(x)
}
\references{
  K. Brendel, E. Comets, C. Laffont, C. Laveille, and
  F.Mentre. Metrics for external model evaluation with an
  application to the population pharmacokinetics of
  gliclazide. \emph{Pharmaceutical Research}, 23:2036--49,
  2006.
}
\seealso{
  \code{\link{set.npde.plotoptions}}
}
\keyword{plot}

