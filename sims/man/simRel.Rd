\name{simRel}
\alias{simRel}
\title{Semantic similarity of Schlicker \emph{et al.} for each pair of terms}
\description{Computes the semantic similarity proposed by Schlicker \emph{et al.} for each pair of terms}
\usage{simRel(sum.mica, ic)}
\arguments{
  \item{sum.mica}{\code{data.frame} generated with function \code{summaryMICA}}
  \item{ic}{\code{numeric} vector with the IC of each term in the ontology}
}
\details{Schlicker \emph{et al.} noted that by taking into account specificities of compared terms can lead to high similarities when comparing general terms. For instance, when comparing general terms the semantic similarity of Lin (see function 'simLin'), the maximal similarity will be obtained comparing a (general) term to itself. Actually, the identity of the indiscernible is generally ensured, with the exception of the root which has an \eqn{IC=0}. However, some treatments require this property not to be respected. Therefore, Schlicker \emph{et al.} introduced a semantic similarity based on a modification of the semantic similarity of Lin, such that

\deqn{sim_{Rel}(t_{i},t_{j}) = sim_{Lin}(t_{i},t_{j})(1 - P(MICA))}
}
\value{The resulting object is a \code{matrix} where for each pair of different terms (rows) is shown the \code{numeric} value of the semantic similarity computed (column)}
\references{
Schlicker, A. \emph{et al.} "A new measure for functional similarity of gene products based on Gene Ontology". BMC Bioinformatics, 7(1):302, 2006.
}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{ancestors}}, \code{\link{resnikSummary}}, \code{\link{simJC}}, \code{\link{simLin}}, \code{\link{simNunivers}}, \code{\link{simPsec}}, \code{\link{simFaith}}, \code{\link{simRes}}, \code{\link{sims.nb}}, \code{\link{summaryMICA}}}
\examples{
data(joslyn)

inv.IminusG <- inverseIminusG(joslyn.OOC)
A.mat <- getA(inv.IminusG)
at <- ancestors(A.mat)
resnik.sum <- resnikSummary(x = joslyn.OOC)
ic <- resnik.sum[, "ic"]
sum.mica <- summaryMICA(at, ic)

sim.Rel <- simRel(sum.mica, ic)
print(sim.Rel)
}
\keyword{ simRel }
\keyword{ ancestors }
\keyword{ getA }
\keyword{ inverseIminusG }
\keyword{ resnikSummary }
\keyword{ summaryMICA }
\keyword{ simRes }
\keyword{ simLin }
\keyword{ simFaith }
\keyword{ simJC }
\keyword{ simNunivers }
\keyword{ simPsec }
\keyword{ sims.nb }
