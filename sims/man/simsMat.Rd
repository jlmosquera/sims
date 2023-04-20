\name{simsMat}
\alias{simsMat}
\title{Coerces a 1-column \code{data.frame} resulting from semantic similarity functions to be an object of class \code{dist}.}
\description{Given a \code{data.frame} with one column resulting from one of the functions for computing semantic similarities, builds a \code{matrix} of class \code{dist}}
\usage{simsMat(x)}
\arguments{
    \item{x}{\code{data.frame} with one column resulting from one of the functions for computing semantic similarities}
}
%\details{}
\value{The resulting object is of class \code{dist}}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{toMat}}, \code{\link{sims.nb}}, \code{\link{sims.eb}}, \code{\link{pseudoDists}}, \code{\link{gosims}}}
\examples{
data(joslyn)

inv.IminusG <- inverseIminusG(joslyn.OOC)
A.mat <- getA(inv.IminusG)
at <- ancestors(A.mat)
resnik.sum <- resnikSummary(x = joslyn.OOC)
ic <- resnik.sum[, "ic"]

sims.Res <- sims.nb(at, ic, method = "Res") 
simsMat(sims.Res)
}
\keyword{ toPairs }
\keyword{ toMat }
\keyword{ sims.nb }
\keyword{ sims.eb }
\keyword{ pseudoDists }
\keyword{ gosims }
\keyword{ ancestors }
\keyword{ resnikSummary }
\keyword{ getA }
\keyword{ inverseIminusG }