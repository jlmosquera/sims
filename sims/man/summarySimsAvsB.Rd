\name{summarySimsAvsB}
\alias{summarySimsAvsB}
\title{Summary of a two-columns \code{matrix} with semantic similarity estimates between each pair of terms for the same measure}
\description{This function provides a summary of estimates calculated with the same semantic similarity measure for two lists of objects (i.e. genes). The summary consists of a \code{list} with three elements: first, a \code{data.frame} with basic statistics for each column, second, the Pearson's Correlation and the correlation test associated with, and three, the Cosine Similarity Measure}
\usage{summarySimsAvsB(x)}
\arguments{
  \item{x}{a \code{matrix} with two columns where for each pair of terms (rows) contains the estimates of the semantic similarity measure for each list of objects (columns)}
}
\details{If any value of the semantic similarities is \code{Inf}, then it is converted to an \code{NA} value, and it is removed from the computation of the statistics.

This function can be used with an arbitrary \code{\link{OOC}} object, but is particularly thought for using Entrex Gene Identifiers (objects) and GO Identifiers (terms).}
%\value{}
%\references{}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{summarySims}}, \code{\link{sims}}, \code{\link{sims.eb}}, \code{\link{pseudoDists}}, \code{\link{gosims}}, \code{\link{gosimsAvsB}}, \code{\link{summarySims}}, \code{\link{cosSim}}}
\examples{
data(prostateIds)                # Data set from the package goProfiles

eg.we <- welsh01EntrezIDs[1:10]   # Entrez Gene IS signature Welsh 01
eg.sg <- singh01EntrezIDs[1:10]   # Entrez Gene IS signature Singh 01

WEvsSG.nb <- gosimsAvsB(eg1 = eg.we, eg2 = eg.sg, ontology = "MF", pckg = "org.Hs.eg.db", type = "nb", method = "Res")

summarySimsAvsB(WEvsSG.nb)
}
\keyword{ summarySimsAvsB }
\keyword{ summarySims }
\keyword{ sims }
\keyword{ sims.eb }
\keyword{ pseudoDists }
\keyword{ gosims }
\keyword{ gosimsAvsB }
\keyword{ cosSim }