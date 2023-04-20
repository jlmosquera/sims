\name{simFaith}
\alias{simFaith}
\title{Semantic similarity of Pirro and Euzenat for each pair of terms}
\description{Computes the semantic similarity proposed by Pirro and Euzenat for each pair of terms}
\usage{simFaith(sum.mica)}
\arguments{
  \item{sum.mica}{\code{data.frame} generated with function \code{summaryMICA}}
}
\details{Pirro and Euzenat introduced a semantic similarity based on the Most Informative Common Ancestor (MICA). It is a node-based approach, such that

\deqn{sim_{Faith}(t_{i},t_{j}) = \frac{IC(MICA)} / {IC(t_{i} + IC(t_{j} - IC(MICA)}}
}
%\value{}
\references{
Pirro, G. and Euzenat, J. "A Feature and Information Theoretic Framework for Semantic Similarity and Relatedness". In Proceedings of the 9th International Semantic Web Conference ISWC. Springer, pp.615-630, 2010.}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{ancestors}}, \code{\link{resnikSummary}}, \code{\link{simJC}}, \code{\link{simLin}}, \code{\link{simNunivers}}, \code{\link{simPsec}}, \code{\link{simRel}}, \code{\link{simRes}}, \code{\link{sims.nb}}, \code{\link{summaryMICA}}}
\examples{
data(joslyn)

inv.IminusG <- inverseIminusG(joslyn.OOC)
A.mat <- getA(inv.IminusG)
at <- ancestors(A.mat)
resnik.sum <- resnikSummary(x = joslyn.OOC)
ic <- resnik.sum[, "ic"]
sum.mica <- summaryMICA(at, ic)

sim.Faith <- simFaith(sum.mica)
print(sim.Faith)
}
\keyword{ simFaith }
\keyword{ ancestors }
\keyword{ getA }
\keyword{ inverseIminusG }
\keyword{ resnikSummary }
\keyword{ summaryMICA }
\keyword{ simRes }
\keyword{ simLin }
\keyword{ simRel }
\keyword{ simJC }
\keyword{ simNunivers }
\keyword{ simPsec }
\keyword{ sims.nb }
