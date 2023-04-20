\name{sims.eb}
\alias{sims.eb}
\title{Wrapper function that calls different methods for computing semantic similarities based on edge-based approaches}
\description{It is wrapper function that calls different methods for computing semantic similarities based on edge-based approaches}
\usage{sims.eb(x, root = NULL, at, method = "Rada")}
\arguments{
  \item{x}{\code{OOC} object}
  \item{root}{\code{charcter} with the name of the root term of the ontology. By default, it takes the first element of 
              vocabulary (i.e. the first element of the slot \code{T} from the \code{OOC} object passed in the argument 
              \code{x}}
  \item{at}{\code{list} of \code{character} vectors with the ancestors of each term}  
  \item{method}{\code{character} indicating the method for computing the semantic similarity. The options implemented in 	
  		\code{sims} package are
               \describe{
                  \item{\code{Rada}}{semantic similarty of Rada}
                  \item{\code{Res.eb}}{semantic similarty of Resnik based on the shortest path (see function \code{\link{distRada}}}
                  \item{\code{all}}{all semantic similarities of edge-based approaches implemented in the package}
                }
             By default computes the semantic similarity of Rada, that is, \code{Rada}}
}
\details{This function computes either all the semantic similarities implemented in the package or the one indicated.}
\value{The resulting object is an \code{data.frame} where for each pair of different terms (rows) are shown the \code{numeric} values of the semantic similarities computed (columns)}
%\references{}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{ancestors}}, \code{\link{simRada}}, \code{\link{simRes.eb}}, \code{\link{sims.nb}}, \code{\link{pseudoDists}}}
\examples{
data(joslyn)

inv.IminusG <- inverseIminusG(joslyn.OOC)
A.mat <- getA(inv.IminusG)
at <- ancestors(A.mat)

sims.eb(x = joslyn.OOC, root = "R", at, method = "Rada")
sims.eb(x = joslyn.OOC, root = "R", at, method = "Res.eb")
sims.eb(x = joslyn.OOC, root = "R", at, method = "all")
}
\keyword{ sims.eb }
\keyword{ ancestors }
\keyword{ summaryMICA }
\keyword{ simRada }
\keyword{ simRes.eb }
\keyword{ sims.nb }
\keyword{ pseudoDists }
