\name{summarySims}
\alias{summarySims}
\title{Summary of semantic similarity estimates between each pair of terms and measure}
\description{For estimates of each measure (columns) provided on matrix, this function builds a \code{data.frame} with the number of pairs, the number of \code{NA}'s, the minimum value and the number of minimum values, the maximum value and the number of minimum values, the mean, the standard deviation and the median}
\usage{summarySims(x)}
\arguments{
  \item{x}{\code{matrix} with the semantic similarity estimates}
}
\details{If any value of the semantic similarities is \code{Inf}, then it is converted to an \code{NA} value, and it is removed from the computation of the statistics}
\value{The resulting object is a \code{data.frame} that for each measure (rows) shows the estimates of statistics provided (columns)}
%\references{}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{ancestors}}, \code{\link{inverseIminusG}}, \code{\link{getA}}, \code{\link{resnikSummary}}, \code{\link{sims}}, \code{\link{sims.eb}}, \code{\link{pseudoDists}}, \code{\link{gosims}}, \code{\link{gosimsAvsB}}, \code{\link{summarySimsAvsB}}}
\examples{
## An arbitrary OOC

data(joslyn)

inv.IminusG <- inverseIminusG(joslyn.OOC)
A.mat <- getA(inv.IminusG)
at <- ancestors(A.mat)
resnik.sum <- resnikSummary(x = joslyn.OOC)
ic <- resnik.sum[, "ic"]
     
sims.all <- sims.nb(at, ic, method = "all")
summarySims(sims.all)

## An OOC associated with the GO

data(prostateIds)                # Data set from the package goProfiles

eg.we <- welsh01EntrezIDs[1:10]   # Entrez Gene IS signature Welsh 01

all.nb <- gosims(eg = eg.we, ontology = "MF", pckg = "org.Hs.eg.db", type = "nb", method = "all")

summarySims(as.matrix(all.nb))
}
\keyword{ summarySims }
\keyword{ inverseIminusG }
\keyword{ getA }
\keyword{ ancestors }
\keyword{ resnikSummary }
\keyword{ sims }
\keyword{ sims.eb }
\keyword{ pseudoDists }
\keyword{ gosims }
\keyword{ gosimsAvsB }
\keyword{ summarySimsAvsB }
