\name{termPairs}
\alias{termPairs}
\title{Builds the pairs of different terms or characters}
\description{Builds a \code{character} vector whose elements are the pairs of different terms or characters}
\usage{termPairs(x)}
\arguments{
	\item{x}{it can be an object of class \code{OOC} or \code{character} vector with the names of each term}
}
\details{Given a \code{character} vector, builds all the pairs of different elements in the vector}
%\value{}
%\references{}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{OOC}}}
\examples{
termPairs(letters[1:5])

data(joslyn)
termPairs(joslyn.OOC)
}
\keyword{ termPairs }
\keyword{ OOC }