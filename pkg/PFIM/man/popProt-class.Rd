\name{popProt-class}
\Rdversion{1.1}
\docType{class}
\alias{popProt-class}
\title{Class \code{"popProt"}}
\description{
  Reference class representing a population protocol
}
\section{Extends}{
All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.
}
%\references{}
%\author{}
%\note{}
\seealso{\code{\link{popProt}}}
%\examples{}
\keyword{classes}
\section{Fields}{
  \describe{
    \item{\code{prots}:}{A list of individual protocols.  Each
      individual protocol consists of a vector of times and a Fisher
      information matrix}
    \item{\code{Ptr}:}{\code{externalptr} to the instance of the C++ class}
  }
}
\section{Methods}{
  \describe{
    \item{\code{Fim()}:}{returns the Fisher information matrix for the protocol}
    \item{\code{LMat()}:}{returns the Cholesky factor of the Fisher
      information matrix}
    \item{\code{ajout()}:}{One additive step in the Fedorov-Wynn
      algorithm.  Returns \code{TRUE} if the step was successful,
      otherwise \code{FALSE}.}
    \item{\code{freq()}:}{returns the vector of frequencies}
    \item{\code{iFim(i, base1)}:}{returns the Fisher information matrix
      for the ith protocol. \code{base1} is a logical scalar indicating
      if the index, \code{i}, is 1-based, as in R.  The alternative is
      0-based, as in C/C++.  By default 1-based indexing is used.}
    \item{\code{iTime(i, base1)}:}{returns the time vector for the ith
      protocol. \code{base1} is a logical scalar indicating
      if the index, \code{i}, is 1-based, as in R.  The alternative is
      0-based, as in C/C++.  By default 1-based indexing is used.}
    \item{\code{ldet()}:}{returns the logarithm of the determinant of
      the Fisher information matrix}
    \item{\code{normalize()}:}{normalize the proportions and drop those
      below the threshold}
    \item{\code{nprot()}:}{returns the number of protocols available}
    \item{\code{ptr()}:}{returns the external pointer, regenerating if
      necessary}
    \item{\code{props()}:}{return the nonzero proportions}
    \item{\code{threshold()}:}{return the threshold for declaring a
      frequency to be zero.} 
    \item{\code{setThreshold(thresh)}:}{set the threshold for declaring
      a proportion to be zero}
  }
}
