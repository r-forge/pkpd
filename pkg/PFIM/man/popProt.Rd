\name{popProt}
\alias{popProt}
\docType{data}
\title{Generator object for the \code{\linkS4class{popProt}} class}
\section{Methods}{
  \describe{
    \item{\code{new(prots=prots)}:}{Create a new \code{\linkS4class{popProt}} object}
  }
}
\arguments{
  \item{prots}{a list of individual protocols.  Each individual protocol
    consist of a vector of times and a Fisher Information matrix, which
    is a positive-definite symmetric matrix.}
}
\note{Arguments to the \code{new} methods must be named arguments.}
\description{The generator object for the \code{\linkS4class{popProt}}
  reference class.  Such an object is primarily used through its
  \code{new} methods.}
\seealso{\code{\linkS4class{popProt}}}
%\examples{}
\keyword{classes}
