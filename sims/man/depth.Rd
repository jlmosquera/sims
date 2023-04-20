\name{depth}
\alias{depth}
\title{Depth of the ontology}
\description{Computes the depth of the ontology. That is, the longest path from the root term to the farthest refinement term end}
\usage{depth(x)}
\arguments{
  \item{x}{can be an \code{OOC} object, a \code{matrix} of the refinements between terms, or a \code{list} of matrices with the number of paths between each pair of terms that are directly connected for each length}
}
%\details{}
\value{\code{numeric} value indicating the longest path from the root term to the far refinement term end.}
%\references{}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{getGk}}, \code{\link{getGr}}, \code{\link{toOOC}}}
\examples{
data(joslyn)

print(mat.g)
d.G <- depth(mat.g)

print(d.G)

## list of matrices with the number of paths
## between each pair of terms for each length k

Gk <- getGk(joslyn.OOC)  
d.Gk <- depth(Gk)

print(d.Gk)

## matrix with all the number of paths of any length
## between each pair of terms

Gr <- getGr(Gk)
d.Gr <- depth(Gr)

print(d.Gr)
}
\keyword{ depth }
\keyword{ getGk }
\keyword{ getGr }
\keyword{ toOOC }
