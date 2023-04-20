\name{summaryPaths}
\alias{summaryPaths}
\title{Lengths of the chains (in terms of depth) or number of paths between each pair of terms.}
\description{Builds a \code{list} of \code{numeric} vectors with either the lengths of the chains (in terms of depth) or the number of paths between each pair of terms.}
\usage{summaryPaths(x, root = NULL, len = TRUE)}
\arguments{
  \item{x}{\code{OOC} object}
  \item{root}{\code{character} indicating the name of the root term of the ontology in the OOC object}
  \item{len}{\code{logical} value. If \code{TRUE}, computes the length (depth) of the chain. If \code{FALSE}, computes the number
  	         of chains. By default is \code{TRUE}}
}
%\details{}
%\value{}
%\references{}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{OOC}}}
\examples{
data(joslyn)

summaryPaths(x = joslyn.OOC, root = "R", len = TRUE)         # Length (depth) of chains
summaryPaths(x = joslyn.OOC, root = "R", len = FALSE)        # Number of chains
}
\keyword{ summaryPaths }
\keyword{ OOC }
