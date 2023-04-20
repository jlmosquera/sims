\name{summaryMICA}
\alias{summaryMICA}
\title{Computes for each pair of terms the Information Content (IC) of each term the Most Informative Common Ancestor (MICA), and the subsumer associated with the MICA}
\description{Builds a \code{data.frame} that for each pair of terms gives the Information Content (IC) of each term the Most Informative Common Ancestor (MICA), and the subsumer associated with the MICA}
\usage{summaryMICA(at, ic)}
\arguments{
  \item{at}{list whose elements are the ancestors of each term}
  \item{ic}{numeric vector with the IC of each term in the ontology}
}
%\details{}
%\value{}
%\references{}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{ancestors}}, \code{\link{commonAncestors}}, \code{\link{inverseIminusG}}, \code{\link{getA}}, \code{\link{resnikSummary}}}
\examples{
data(joslyn)

inv.IminusG <- inverseIminusG(joslyn.OOC)
A.mat <- getA(inv.IminusG)
at <- ancestors(A.mat)
ca <- commonAncestors(at)
resnik.sum <- resnikSummary(x = joslyn.OOC)
ic <- resnik.sum[, "ic"]

sum.mica <- summaryMICA(at, ic)
print(sum.mica)
}
\keyword{ summaryMICA }
\keyword{ ancestors }
\keyword{ inverseIminusG }
\keyword{ getA }
\keyword{ commonAncestors }
\keyword{ resnikSummary }
