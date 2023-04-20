\name{ICA}
\alias{ICA}
\title{Information Content (IC) of common ancestors}
\description{Builds a \code{list} of numeric vectors with the Information Content (IC) associated with each common ancestor of each pair of terms}
\usage{ICA(x, ic)}
\arguments{
  \item{x}{\code{list} of common ancestors for each pair of terms}
  \item{ic}{\code{numeric} vector with the IC of each term in the ontology}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The resulting object is a \code{list} whose elements are \code{numeric} vectors with the ICs of the common ancestors for each pair of terms}
\references{
Resnik P. "Using information content to evaluate semantic similarity in a taxonomy". Proceedings of the 14th International Joint Conference on Artificial Intelligence, pp. 448-453, 1995.}
\author{Jose Luis Mosquera}
%\note{
%%  ~~further notes~~
%}
\seealso{\code{\link{ancestors}}, \code{\link{commonAncestors}}, \code{\link{inverseIminusG}}, \code{\link{resnikSummary}}
}
\examples{
data(joslyn)

inv.IminusG <- inverseIminusG(joslyn.OOC)     # (I - G)^(-1)
A.mat <- getA(inv.IminusG)                    # accessibility matrix
at <- ancestors(A.mat)                        # ancestors
ca <- commonAncestors(at)                     # common ancestors
resnik.sum <- resnikSummary(x = joslyn.OOC)   # summary of Resnik measures

ica <- ICA(x = ca, ic = resnik.sum[, "ic"])   
print(ica)
}
\keyword{ ICA }
\keyword{ ancestors }
\keyword{ commonAncestors }
\keyword{ resnikSummary }
