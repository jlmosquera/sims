\name{simRada}
\alias{simRada}
\title{Semantic similarity measure of Rada \emph{et al.} for each pair of terms}
\description{Computes the semantic similarity, based on the shortest path, proposed by Rada \emph{et al.} between each pair of terms}
\usage{simRada(sum.paths, at)}
\arguments{
  \item{sum.paths}{\code{list} of \code{numeric} vectors with the lengths (in terms of depth) of the number of paths 
                   between each pair of terms.}
  \item{at}{\code{list} of ancestors of each term}
}
\details{The measure involved in this function is an edge-based approach. It was proposed by Rada \emph{et al.}. The measure is a transformation from the distance of Rada \emph{et al.} (see function \code{\link{distRada}}) to semantic similarity through the formula \eqn{sim = \frac{1}{1 + dist}}. That is, the semantic similarity proposed by Rada \emph{et al.} is defined as

\deqn{sim_{Rada}(t_{i}, t_{j})=\frac{1}{1 + d_{Rada}(t_{i}, t_{j})}}
}
\value{The resulting object is a \code{matrix} where for each pair of different terms (rows) is shown the \code{numeric} value of the semantic similarity computed (column).}
\references{
Rada, R. et al. Development and application of a metric on semantic nets. Ieee Transactions On Systems Man And Cybernetics, 19(1), pp.17-30, 1989.
}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{ancestors}}, \code{\link{simRes.eb}}, \code{\link{distRada}}, \code{\link{LCAs}}, \code{\link{summaryPaths}}}
\examples{
data(joslyn)

sum.paths <- summaryPaths(x = joslyn.OOC, root = "R", len = TRUE)

inv.IminusG <- inverseIminusG(joslyn.OOC)
A.mat <- getA(inv.IminusG)
at <- ancestors(A.mat)

sim.Rada <- simRada(sum.paths, at)
print(sim.Rada)
}
\keyword{ simRada }
\keyword{ simRes.eb }
\keyword{ distRada }
\keyword{ ancestors }
\keyword{ LCAs }
\keyword{ summaryPaths }
\keyword{ getA }
\keyword{ inverseIminusG }