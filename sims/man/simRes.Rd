\name{simRes}
\alias{simRes}
\title{Semantic similarity of Resnik for each pair of terms}
\description{Computes the semantic similarity proposed by Resnik for each pair of terms}
\usage{simRes(at, ic, subs = FALSE)}
\arguments{
  \item{at}{ list whose elements are the ancestors of each term}
  \item{ic}{numeric vector with the IC of each term in the ontology}
  \item{subs}{logical value. If TRUE, then the resulting data.frame will show an etra column with the name of the subsumer
              term associated with MICA value}
}
\details{Given two terms of an ontology the Most Informative Common Ancestor (MICA) is the common ancestor with the highest Information Content (IC). This measure is a semantic similarity, from the node-based approaches, introduced by Resnik and defined as

\deqn{sim_{Res}(t_{i},t_{j}) = IC(MICA) = \max_{t\in S(t_{i},t{j})} (IC(t)}

where \eqn{S(t_{i,},t{j})} is the set of terms that subsumes both terms \eqn{t_{i}} and \eqn{t_{j}}, and \eqn{IC(t)} is the information content measure of \eqn{t}.

MICA does not take into account the disjoint common ancestors, that is, those common ancestors that do not subsume any other common ancestor.}
\value{The resulting object is a \code{matrix} where for each pair of different terms (rows) is shown the \code{numeric} value of the semantic similarity computed (column)}
\references{
Resnik P. "Using information content to evaluate semantic similarity in a taxonomy". Proceedings of the 14th International Joint Conference on Artificial Intelligence. 448-453, 1995
}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{ancestors}}, \code{\link{commonAncestors}}, \code{\link{resnikSummary}}, \code{\link{simJC}}, \code{\link{simLin}}, \code{\link{simNunivers}}, \code{\link{simPsec}}, \code{\link{simRel}}, \code{\link{simFaith}}, \code{\link{sims.nb}}, \code{\link{summaryMICA}}}
\examples{
data(joslyn)

inv.IminusG <- inverseIminusG(joslyn.OOC)
A.mat <- getA(inv.IminusG)
at <- ancestors(A.mat)
ca <- commonAncestors(at)
resnik.sum <- resnikSummary(x = joslyn.OOC)

sim.Res <- simRes(at, ic = resnik.sum[, "ic"])
print(sim.Res)
}
\keyword{ simRes }
\keyword{ ancestors }
\keyword{ commonAncestors }
\keyword{ getA }
\keyword{ inverseIminusG }
\keyword{ resnikSummary }
\keyword{ summaryMICA }
\keyword{ simFaith }
\keyword{ simLin }
\keyword{ simRel }
\keyword{ simJC }
\keyword{ simNunivers }
\keyword{ simPsec }
\keyword{ sims.nb }
