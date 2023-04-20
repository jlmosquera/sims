\name{getGr}
\alias{getGr}
\title{Number of paths of any length between each pair of terms that are directly connected}
\description{Computes all the number of paths of any length between each pair of terms that are directly connected.}
\usage{getGr(x)}
\arguments{
  \item{x}{can be an object \code{OOC}, a \code{matrix} of refinement matrix, or a \code{list} of matrices with the number of 
           paths of each length between each pair of terms that are directly connected}
}
\details{The computation of the number of paths is based on the formula 

\deqn{\mathbf{I} + \Gamma^1 + \ldots + \Gamma^r = (\mathbf{I} - \Gamma)^(-1)}.

There are two ways for performing the computation; first, by calculating the sum of the power matrices of the refinements \code{matrix}, and second, by computing the inverse of the \code{matrix} \eqn{\mathbf{I} - \Gamma}.

Note that depending on the class of the argument \code{x} supplied, the computation requires calculating the power matrices of the refinements matrix. That is, \eqn{\Gamma^1, \Gamma^2,\ldots, \Gamma^r}, where \code{r} is the depth of the ontology}
\value{Resulting object is a \code{matrix} whose elements are the number of paths of any length between each pair of terms that are directly connected}
%\references{}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{getGk}}, \code{\link{toOOC}}}
\examples{
data(joslyn)

## Based on the OOC object

Gr <- getGr(joslyn.OOC)
print(Gr)

## Based on the refinement matrix

Gr.2 <- getGr(mat.g)
Gr==Gr.2

## Based on the list of matrices with the number of 
## paths of each length between each pair of terms 
## that are directly connected

Gk <- getGk(joslyn.OOC)
Gr.3 <- getGr(Gk)
Gr==Gr.3
}
\keyword{ getGr }
\keyword{ getGk }
\keyword{ toOOC }
