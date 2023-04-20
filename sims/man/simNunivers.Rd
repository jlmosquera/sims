\name{simNunivers}
\alias{simNunivers}
\title{Semantic similarity of Mazandu and Mulder for each pair of terms}
\description{Computes the semantic similarity proposed by Mazandu and Mulder for each pair of terms}
\usage{simNunivers(sum.mica)}
\arguments{
  \item{sum.mica}{\code{data.frame} generated with function \code{summaryMICA}}
}
\details{Mazandu and Mulder introduced a semantic similarity based on the Most Informative Common Ancestor (MICA). It is a node-based approach, such that

\deqn{sim_{Nunivers}(t_{i},t_{j}) = \frac{IC(MICA)}{\max \{IC(t_{i}), IC(t_{j}\}}}
}
%\value{}
\references{
Mazandu, G.K. and Mulder, N.J. "Information content-based Gene Ontology semantic similarity approaches: Toward a unified framework theory". BioMed Research International, 2013.
}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{ancestors}}, \code{\link{resnikSummary}}, \code{\link{simJC}}, \code{\link{simLin}}, \code{\link{simFaith}}, \code{\link{simPsec}}, \code{\link{simRel}}, \code{\link{simRes}}, \code{\link{sims.nb}}, \code{\link{summaryMICA}}}
\examples{
data(joslyn)

inv.IminusG <- inverseIminusG(joslyn.OOC)
A.mat <- getA(inv.IminusG)
at <- ancestors(A.mat)
resnik.sum <- resnikSummary(x = joslyn.OOC)
ic <- resnik.sum[, "ic"]
sum.mica <- summaryMICA(at, ic)

sim.Nunivers <- simNunivers(sum.mica)
print(sim.Nunivers)
}
\keyword{ simNunivers }
\keyword{ ancestors }
\keyword{ getA }
\keyword{ inverseIminusG }
\keyword{ resnikSummary }
\keyword{ summaryMICA }
\keyword{ simRes }
\keyword{ simLin }
\keyword{ simRel }
\keyword{ simJC }
\keyword{ simFaith }
\keyword{ simPsec }
\keyword{ sims.nb }
