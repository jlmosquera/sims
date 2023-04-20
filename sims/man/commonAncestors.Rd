\name{commonAncestors}
\alias{commonAncestors}
\title{Common ancestors for each pair of terms of the ontology}
\description{Generates a list whose elements are the common ancestors for each pair of terms of the ontology}
\usage{commonAncestors(at)}
\arguments{
	\item{at}{\code{list} of \code{character} vectors containing the name of the ancestors of each term of
	          the ontology}
}
%\details{}
\value{Resulting object is a \code{list} whose elements are \code{character} vectors with the names of the common	ancestors between each pair of terms connected}
%\references{}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{ancestors}}, \code{\link{getA}}, \code{\link{inverseIminusG}}
}
\examples{
data(joslyn)

inv.IminusG <- inverseIminusG(joslyn.OOC)
A.mat <- getA(inv.IminusG)
at <- ancestors(A.mat)

ca <- commonAncestors(at)
print(at)
}
\keyword{ commonAncestors }
\keyword{ inverseIminusG }
\keyword{ geatA }
\keyword{ ancestors }
