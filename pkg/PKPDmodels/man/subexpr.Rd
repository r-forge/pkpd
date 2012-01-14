\name{subexpr}
\alias{subexpr}
\title{"Sub[stitute] expr[ession]"}
\usage{
  subexpr(expr, nm, sub)
}
\arguments{
  \item{expr}{an expression}

  \item{nm}{a name for which to substitute}

  \item{sub}{the expression to substitute for name nm}
}
\value{
  the expression with all occurrences of nm replaced by sub
}
\description{
  Substitute the expression \code{sub} for the name
  \code{nm} in \code{expr} by walking the tree.
}
\note{
  this function is called recursively
}

