\name{simJC}
\alias{simJC}
\title{Semantic similarity of Jiang and Conrath for each pair of terms}
\description{Computes the semantic similarity proposed by Jiang and Conrath for each pair of terms}
\usage{simJC(sum.mica)}
\arguments{
  \item{sum.mica}{\code{data.frame} generated with function \code{summaryMICA}}
}
\details{Jiang and Conrath introduced a semantic similarity based on the Most Informative Common Ancestor (MICA). It is a node-based approach, such that

\deqn{sim_{JC}(t_{i},t_{j}) = \frac{1}{1 + dist_{JC}(t_{i},t_{j})} = \frac{1}{1 + (IC(t_{i}) + IC(t_{j}) - (2 MICA))}}
}
%\value{}
\references{
Jiang, J.J. and Conrath, D.W. "Semantic similarity based on corpus statistics and lexical taxonomy". In Proceedings of International Conference Research on Computational Linguistics (ROCLING X), pp.19-33, 1997.
}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{ancestors}}, \code{\link{resnikSummary}}, \code{\link{simFaith}}, \code{\link{simLin}}, \code{\link{simNunivers}}, \code{\link{simPsec}}, \code{\link{simRel}}, \code{\link{simRes}}, \code{\link{sims.nb}}, \code{\link{summaryMICA}}}
\examples{
data(joslyn)

inv.IminusG <- inverseIminusG(joslyn.OOC)
A.mat <- getA(inv.IminusG)
at <- ancestors(A.mat)
resnik.sum <- resnikSummary(x = joslyn.OOC)
ic <- resnik.sum[, "ic"]
sum.mica <- summaryMICA(at, ic)

sim.JC <- simJC(sum.mica)
print(sim.JC)
}
\keyword{ simJC }
\keyword{ ancestors }
\keyword{ getA }
\keyword{ inverseIminusG }
\keyword{ resnikSummary }
\keyword{ summaryMICA }
\keyword{ simRes }
\keyword{ simLin }
\keyword{ simRel }
\keyword{ simFaith }
\keyword{ simNunivers }
\keyword{ simPsec }
\keyword{ sims.nb }
