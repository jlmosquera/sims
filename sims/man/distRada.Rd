\name{distRada}
\alias{distRada}
\title{Distances of the shortest paths between each pair of terms in the ontology}
\description{This function is a edge-based approach to compute the distance of the shortest between each pair of terms in the ontology.}
\usage{distRada(sum.paths, at)}
\arguments{
  \item{sum.paths}{\code{list} of \code{numeric} vectors with the lengths (in terms of depth) of the number of paths 
                   between each pair of terms.}
  \item{at}{\code{list} of ancestors of each term}
}
\details{This function is based on the edge-counting method (i.e. edge-based approach) proposed by Rada \emph{et al.} The measure computes the length of the shortest path (i.e. the number of edges between a pair of terms). 

This taxonomical distance was thought to be applied in trees where the shortest path between two terms contains a unique common ancestor, which is known as the Least Common Ancestor (LCA). Note that, in ontologies, this measure is not the most appropriate to be used because DAG structures may show more than one LCA. In any case, the method implemented does not distinguish a unique LCA (see function \code{\link{LCAs}})}
\value{Resulting object is a \code{matrix} whose elements show the distance (column) between each pair of different terms (rows)}
\references{
Rada, R. et al. Development and application of a metric on semantic nets. Ieee Transactions On Systems Man And Cybernetics, 19(1), pp.17-30, 1989.
}
\author{Jose Luis Mosquera}
%\note{}

\seealso{\code{\link{ancestors}}, \code{\link{LCAs}}, \code{\link{simRada}},  \code{\link{simRes.eb}}, \code{\link{summaryPaths}}}
\examples{
data(joslyn)

sum.paths <- summaryPaths(x = joslyn.OOC, root = "R", len = TRUE)

inv.IminusG <- inverseIminusG(joslyn.OOC)
A.mat <- getA(inv.IminusG)
at <- ancestors(A.mat)

dist.rada <- distRada(sum.paths,at)
print(dist.rada)
}
\keyword{ distRada }
\keyword{ ancestors }
\keyword{ inverseIminusG }
\keyword{ getA }
\keyword{ LCAs }
\keyword{ simRada }
\keyword{ simRes.eb }
\keyword{ summaryPaths }
\keyword{ getA }
\keyword{ inverseIminusG }
