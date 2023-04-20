\name{simRes.eb}
\alias{simRes.eb}
\title{Semantic similarity measure of Resnik \emph{et al.} for each pair of terms, considering the maximal depth of the ontology}
\description{Computes the semantic similarity measure of Resnik between each pair of node, considering the maximal depth of the ontology.}
\usage{simRes.eb(sum.paths, at, x)}
\arguments{
  \item{sum.paths}{\code{list} of \code{numeric} vectors with the lengths (in terms of depth) of the number of paths 
                   between each pair of terms.}
  \item{at}{\code{list} of ancestors of each term}
  \item{x}{\code{OOC} object}
}
\details{This function computes a measure proposed by Resnik in order to normalize the measure of semantic evidences from the graph of the ontology and computed with the the shortest path containing the LCAs of the terms involved int (see function \code{\link{distRada}}). It is an edge-based approach, and given two terms \eqn{t_{i}} and \eqn{t_{j}}, it is defined as  

\deqn{latex}{sim_{Res-ed}(t_{i},t_{j}) = 2 * max.depth - (dist_{Rada}(t_{i},t_{LCA}) + dist_{Rada}(t_{j},t_{LCA}))}
}
\value{The resulting object is a \code{matrix} where for each pair of different terms (rows) is shown the \code{numeric} value of the semantic similarity computed (column).}
\references{
Resnik P. "Using information content to evaluate semantic similarity in a taxonomy". Proceedings of the 14th International Joint Conference on Artificial Intelligence. 448-453, 1995

Rada, R. et al. Development and application of a metric on semantic nets. Ieee Transactions On Systems Man And Cybernetics, 19(1), pp.17-30, 1989.
}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{distRada}}, \code{\link{simRada}}}
\examples{
data(joslyn)

sum.paths <- summaryPaths(x = joslyn.OOC, root = "R", len = TRUE)

inv.IminusG <- inverseIminusG(joslyn.OOC)
A.mat <- getA(inv.IminusG)
at <- ancestors(A.mat)

sim.Res.eb <- simRes.eb(sum.paths, at, x = joslyn.OOC)
print(sim.Res.eb)
}
\keyword{ simRes.eb }
\keyword{ simRada }
\keyword{ distRada }
