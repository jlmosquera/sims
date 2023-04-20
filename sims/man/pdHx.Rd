\name{pdHx}
\alias{pdHx}
\title{Pseudo-distance of the maximum chain lengths between comparable terms of the ontology.}
\description{This function computes pseudo-distances of the maximum chain lengths proposed by Josly \emph{et al.}}
\usage{pdHx(sum.paths)}
\arguments{
  \item{sum.paths}{\code{list} of vectors with the lengths of the chains of each pair of terms}
}
\details{Pseudo-distance concept is an strategy of the edge-based approach for measuring the distances between terms in an ontology. It was proposed by Joslyn \emph{et al.}

\code{pdHx} computes the pseudo-distances proposed of the maximum chain lengths. That is, let \eqn{t_{i}} and \eqn{t_{j}} be two comparable terms, then the maximum chain lengths is defined as

\deqn{\delta_{x}(t_{i},t_{j})=h^{*}(t_{i},t_{j})=\max_{C\in\mathcal{C}(t_{i},t_{j})}C\vert C\vert}

In case of non-comparable terms, the resulting value is NA.}
%\value{}
\references{
Cliff A. Joslyn, Susan M. Mniszewski Andy W. Fulmer, and Gary G. Heaton. (2004). "The gene ontology categorizer". Bioinformatics, 20(s1):169-77, 2004.}
\author{Jose Luis Mosquera}
%\note{
%%  ~~further notes~~
%}
\seealso{\code{\link{pdHm}}, \code{\link{pdHax}}, \code{\link{pdHap}}, \code{\link{pseudoDists}}}
\examples{
data(joslyn)

length.chains <- summaryPaths(x = joslyn.OOC, root = "R", len = TRUE)   # length of chains
num.chains <- summaryPaths(x = joslyn.OOC, root = "R", len = FALSE)      # number of chains                ##                                                   
sum.paths <- lapply(1:nrow(length.chains), function(i, x, y)
                                           {
                                               ch <- rep(x[i, ], y[i, ])
                                               if(sum(ch)==0) ch <- 0
                                               return(ch)
                                           },
                                           x = length.chains,
                                           y = num.chains)
names(sum.paths) <- rownames(length.chains)

pseudodist.Hx <- pdHx(sum.paths)
print(pseudodist.Hx)
}
\keyword{ pdHx }
\keyword{ pdHm }
\keyword{ pdHax }
\keyword{ pHap }
\keyword{ pseudoDists }
