\name{getA}
\alias{getA}
\title{Accessibility matrix associated with the DAG structure of the ontology}
\description{Builds the accessibility matrix associated with the DAG structure of the ontology}
\usage{getA(x)}
\arguments{
	\item{x}{\code{matrix} of the number of paths of any length between each pair of terms.}
}
%\details{The argument \code{x} accepts two possible matrices with the number of paths of any between each pair
%of terms. On one hand, can be supplied the resulting \code{matrix} of the function \code{inverseIminusG}, and
%on the other hand, the resulting \code{matrix} generated with the function \code{getGr}. }

\value{Resulting object is a \code{matrix} whose elements contains a \code{TRUE} value, if there exists a path of terms path following the arcs from the term on the rows to the term on the columns.

Depending on the approach applied later on, resulting \code{matrix} may be interpreted as the accessibility matrix (i.e. a node-based approach, and an edge-based approach based on distances, or a hybrid-based approach) or the comparability matrix (i.e. an edge-based based on pseudo-distances). In the first case, \code{TRUE} values indicate terms that are connected by a path, and \code{FALSE} values indicate which not. In the second case, \code{TRUE} values indicate comparable nodes (i.e. terms),
but \code{FALSE} values indicate non-comparable nodes.}
%\references{}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{inverseIminusG}}, \code{\link{getGr}}}
\examples{
data(joslyn)

## (I-Gamma)^{-1}
inv.IminusG <- inverseIminusG(joslyn.OOC)
getA(inv.IminusG)

## I + Gamma^1 + Gamma^2 + ... + Gamma^r, where r is the depth of the ontology.
Gr <- getGr(joslyn.OOC)

## (I-Gamma)^{-1} = I + Gamma^1 + Gamma^2 + ... + Gamma^r                        
Gr==inv.IminusG-diag(nrow(inv.IminusG))        

## Accessibility matrix
getA(inv.IminusG-diag(nrow(inv.IminusG)))

getA(Gr)
}
\keyword{ getA }
\keyword{ inverseIminusG }
\keyword{ getGr }
