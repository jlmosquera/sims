\name{simLin}
\alias{simLin}
\title{Semantic similarity of Lin for each pair of terms}
\description{Computes the semantic similarity proposed by Lin for each pair of terms}
\usage{simLin(sum.mica)}
\arguments{
  \item{sum.mica}{\code{data.frame} generated with function \code{summaryMICA}}
}
\details{Lin introduced a semantic similarity based on the Most Informative Common Ancestor (MICA). It is a node-based approach, such that

\deqn{sim_{Lin}(t_{i},t_{j}) = 2 IC(MICA) / (IC(t_{i}) + IC(t_{j}))}
}
%\value{}
\references{
Lin, D. "An information-theoretic definition of similarity". In Proceedings of the Fifteenth International Conference on Machine Learning, Morgan Kaufmann Publishers, pp. 296–304, 1998.
}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{ancestors}}, \code{\link{resnikSummary}}, \code{\link{simJC}}, \code{\link{simFaith}}, \code{\link{simNunivers}}, \code{\link{simPsec}}, \code{\link{simRel}}, \code{\link{simRes}}, \code{\link{sims.nb}}, \code{\link{summaryMICA}}}
\examples{
data(joslyn)

inv.IminusG <- inverseIminusG(joslyn.OOC)
A.mat <- getA(inv.IminusG)
at <- ancestors(A.mat)
resnik.sum <- resnikSummary(x = joslyn.OOC)
ic <- resnik.sum[, "ic"]
sum.mica <- summaryMICA(at, ic)

sim.Lin <- simLin(sum.mica)
print(sim.Lin)
}
\keyword{ simLin }
\keyword{ ancestors }
\keyword{ getA }
\keyword{ inverseIminusG }
\keyword{ resnikSummary }
\keyword{ summaryMICA }
\keyword{ simRes }
\keyword{ simFaith }
\keyword{ simRel }
\keyword{ simJC }
\keyword{ simNunivers }
\keyword{ simPsec }
\keyword{ sims.nb }
