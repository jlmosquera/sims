\name{LCAs}
\alias{LCAs}
\title{Length of the shortest paths containing the Least Common Ancestors (LCA) between each pair of terms}
\description{Computes the Least Common Ancestors between each pair of nodes}
\usage{LCAs(at, sum.paths)}
\arguments{
  \item{at}{\code{list} of ancestors of each term}
  \item{sum.paths}{\code{list} of \code{numeric} vectors with the lengths (in terms of depth) of the number of paths 
                   between each pair of terms.}
}
\details{The function looks for all the least common ancestor of each pair of terms and computes the shortest path between such terms containing the corresponding LCA.}
\value{The resulting object is a \code{list} whose elements are \code{numeric} vectors with the lengths of the shortest paths associated with each LCA.}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{jose Luis Mosquera}
%\note{
%%  ~~further notes~~
%}
\seealso{\code{\link{ancestors}}, \code{\link{summaryPaths}}
}
\examples{
data(joslyn)

sum.paths <- summaryPaths(x = joslyn.OOC, root = "R", len = TRUE)

inv.IminusG <- inverseIminusG(joslyn.OOC)
A.mat <- getA(inv.IminusG)
at <- ancestors(A.mat)

lcas <- LCAs(at, sum.paths)
print(lcas)
}
\keyword{ LCAs }
\keyword{ ancestors }
\keyword{ inverseIminusG }
\keyword{ getA }
\keyword{ summaryPaths }