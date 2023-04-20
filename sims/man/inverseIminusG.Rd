\name{inverseIminusG}
\alias{inverseIminusG}
\title{Computes the number of paths of any length between each pair of terms in the ontology}
\description{Builds a \code{matrix} with the number of paths of any length between each pair of terms in the ontology}
\usage{inverseIminusG(x)}
\arguments{
  \item{x}{can be an object of class \code{OOC} or a \code{matrix} with the refinements of each term in the ontology}
}
\details{In order to compute faster the matrix, applies the formula \code{(I - Gamma)^(-1)}, instead of using \code{I + Gamma^1 + \ldots + Gamma^r} where \code{r} is the depth of the ontology.}
\value{The resulting object is a \code{matrix} whose elements are the number of paths}
%\references{}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{toOOC}}}
\examples{
data(joslyn)

inv.IminusG <- inverseIminusG(joslyn.OOC)
print(inv.IminusG)
}
\keyword{ inverseIminusG }
\keyword{ toOOC }