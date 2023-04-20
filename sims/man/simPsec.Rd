\name{simPsec}
\alias{simPsec}
\title{Semantic similarity of Pirro and Seco for each pair of terms}
\description{Computes the semantic similarity proposed by Pirro and Seco for each pair of terms}
\usage{simPsec(sum.mica)}
\arguments{
  \item{sum.mica}{\code{data.frame} generated with function \code{summaryMICA}}
}
\details{Pirro and Seco introduced a semantic similarity based on the Most Informative Common Ancestor (MICA). It is a node-based approach, such that

\deqn{sim_{Psec}(t_{i},t_{j}) = (3 IC(MICA)) - IC(t_{i} - IC(t_{j}}
}
%\value{}
\references{
Pirro, G. and Seco, N. "Design, Implementation and Evaluation of a New Semantic Similarity Metric Combining Features and Intrinsic Information Content". Meersman, R. and Tari, Z. eds. Lecture Notes in Computer Science, 5332:1271-1288, 2008.

Pirro, G. "A semantic similarity metric combining features and intrinsic information content". Data and Knowledge Engineering, 68(11):1289-1308, 2009.
}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{ancestors}}, \code{\link{resnikSummary}}, \code{\link{simJC}}, \code{\link{simLin}}, \code{\link{simNunivers}}, \code{\link{simFaith}}, \code{\link{simRel}}, \code{\link{simRes}}, \code{\link{sims.nb}}, \code{\link{summaryMICA}}}
\examples{
data(joslyn)

inv.IminusG <- inverseIminusG(joslyn.OOC)
A.mat <- getA(inv.IminusG)
at <- ancestors(A.mat)
resnik.sum <- resnikSummary(x = joslyn.OOC)
ic <- resnik.sum[, "ic"]
sum.mica <- summaryMICA(at, ic)

sim.Psec <- simPsec(sum.mica)
print(sim.Psec)
}
\keyword{ simPsec }
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
\keyword{ simFaith }
\keyword{ sims.nb }
