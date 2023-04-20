\name{resnikSummary}
\alias{resnikSummary}
\title{Summary table providing with the number of times that each term or any of its refinements appears in the OOC, the probability of finding the term, and the Information Content of the term}
\description{Builds a \code{data.frame} providing the number of times that each term or any of its refinements appears in the OOC (i.e. \code{n(t_{i})}), the probability of finding the term (i.e. \code{p(t_{i})}), and the Information Content of the term (i.e. \code{IC(t_{i})})}
\usage{resnikSummary(x, root = NULL)}
\arguments{
  \item{x}{an object of class \code{OOC}}
  \item{root}{a \code{character} with the name of the root term in the ontology. If \code{NULL}, it takes the first element of 
              vocabulary, that is, the first element of the slot \code{T} from the object \code{OOC} \code{x})}
}
%\details{}
\value{
The resulting object is a \code{data.frame} where on rows are all the terms of the ontology, and on columns are the following measures
\item{nt }{the number of times that each term or any of its refinements appears in the OOC}
\item{pt }{the probability of finding the term}
\item{ic }{the Information Content of the term}
}
\references{
Resnik, P. "Using information content to evaluate semantic similarity in a taxonomy". Proceeding of the 14th International Joint Conference on Artificial Intelligence. pp. 448-453, 1995.
}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{Nt}}, \code{\link{toOOC}}}
\examples{
data(joslyn)

sum.resnik <- resnikSummary(joslyn.OOC)
print(sum.resnik)
}
\keyword{ resnikSummary }
\keyword{ Nt }
\keyword{ toOOC }
