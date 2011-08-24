\name{subexpr}
\alias{subexpr}
\title{"Sub[stitute] expr[ession]"}
\usage{subexpr(expr, nm, sub)}
\description{Substitute the expression \code{sub} for the name \code{nm} in
\code{expr} by walking the tree.}
\value{the expression with all occurrences of nm replaced by sub}
\note{this function is called recursively}
\arguments{\item{expr}{an expression}
\item{nm}{a name for which to substitute}
\item{sub}{the expression to substitute for name nm}}
