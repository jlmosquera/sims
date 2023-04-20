\name{pseudoDists}
\alias{pseudoDists}
\title{Wrapper function that calls different methods for computing pseudo-distances}
\description{It is a wrapper function that calls different methods for computing pseudo-distances}
\usage{pseudoDists(x, root = NULL, method = "hm")}
\arguments{
  \item{x}{\code{OOC} object}
  \item{root}{\code{character} with the name of the root term of the ontology in the \code{OOC} object}
  \item{method}{\code{character} indicating the method for computing the pseudo-distance. The options implemented in \code{sims} package are
               \itemize{
                  \item{hm}{pseudo-distance of the minimum chain length}
                  \item{hx}{pseudo-distance of the maximum chain length}
                  \item{hax}{pseudo-distance of the average of extreme chain lengths}
                  \item{hap}{pseudo-distance of the average of all chain lengths}
                  \item{all}{all pseudo-distances}
               }
               By default computes the pseudo-distance of the minimum chain length, i.e. \code{hm}}
}
\details{In case of non-comparable terms, the resulting value is NA}
\value{The resulting object is a \code{data.frame} where for each pair of different terms (rows) are shown the pseudo-distances proposed by Joslyn \emph{et al.} (columns).}
\references{
Cliff A. Joslyn, Susan M. Mniszewski Andy W. Fulmer, and Gary G. Heaton. (2004). "The gene ontology categorizer". Bioinformatics, 20(s1):169-77, 2004.}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{toOOC}}, \code{\link{pdHm}}, \code{\link{pdHx}}, \code{\link{pdHax}}, \code{\link{pdHap}}}
\examples{
data(joslyn)

pseudoDists(x = joslyn.OOC, root = "R", method = "hm")          # minimum chain length
pseudoDists(x = joslyn.OOC, root = "R", method = "hx")          # maximum chain length
pseudoDists(x = joslyn.OOC, root = "R", method = "hax")         # average of extreme chain lengths
pseudoDists(x = joslyn.OOC, root = "R", method = "hap")         # average of all chain lengths

pseudoDists(x = joslyn.OOC, root = "R", method = "all")         # all pseudo-distances
}
\keyword{ pseudoDists }
\keyword{ pdHm }
\keyword{ pdHx }
\keyword{ pdHax }
\keyword{ pdHap }
\keyword{ toOOC }
