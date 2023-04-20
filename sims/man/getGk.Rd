\name{getGk}
\alias{getGk}
\title{Builds a \code{list} of the matrices with the number of paths between each pair of terms that are directly connected for each length}
\description{Builds a \code{list} whose elements are the matrices (one per each length) of the number of paths between each pair of terms that are directly connected. That is, \eqn{\Gamma^1, \Gamma^2, \ldots, \Gamma^r}, where \code{r} is the depth of the ontology}
\usage{getGk(x)}
\arguments{
  \item{x}{an \code{OOC} or a \code{matrix} of refinements}
}
\details{Given either an \code{\link{OOC}} object or a refinement \code{matrix} (i.e. \code{G} matrix in the \code{\link{OOC}} object), computes all the number of paths for each length between each pair of terms that are directly connected. Each of these matrices are power matrices of the matrix \code{G}. For example, to obtain all the number of paths between each pair of terms connected in the ontology with length 1, the function computes \eqn{D^1}, to obtain all the number of paths with length 2, calculates \eqn{G^2 = G * G}, and so on until the length of the largest path, that is, the maximum depth of the ontology}
\value{The resulting object is a \code{list} of matrices. Each \code{matrix} is associated with a length of paths between each pair of connected terms. Elements of these matrices are the number of paths for the same length connecting each pair of terms. For example, let \code{r} be the depth of the ontology, then the elements of the resulting \code{list} are
\describe{
	\item{\code{G1}}{\code{matrix} with the number of paths between each pair of terms with length 1}
	\item{\code{G2}}{\code{matrix} with the number of paths between each pair of terms with length 2}
	\item{\ldots}{\ldots}
	\item{\code{Gr}}{\code{matrix} with the number of paths between each pair of terms with length \code{r}}
}
}
%\references{}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{getGr}}, \code{\link{depth}}, \code{\link{toOOC}}}
\examples{
data(joslyn)

Gk <- getGk(joslyn.OOC)  # based on an OOC object
print(Gk)

Gk.2 <- getGk(mat.g)         # based on a matrix
print(Gk.2)
}
\keyword{ getGk }
\keyword{ getGr }
\keyword{ depth }
\keyword{ toOOC }
