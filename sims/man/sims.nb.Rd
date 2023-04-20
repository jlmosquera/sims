\name{sims.nb}
\alias{sims.nb}
\title{Wrapper function that calls different methods for computing semantic similarities based on node-based approaches}
\description{It is a wrapper function that calls different methods for computing semantic similarities based on node-based approaches}
\usage{sims.nb(at, ic, method = "Res")}
\arguments{
  \item{at}{\code{list} of \code{character} vectors with the ancestors of each term}
  \item{ic}{\code{numeric} vector with the IC of each term in the ontology}
  \item{method}{\code{character} indicating the method for computing the semantic similarity.The options implemented in 	
                \code{sims} package are
               \describe{
                  \item{\code{Res}}{semantic similarty of Resnik}
                  \item{\code{Lin}}{semantic similarty of Lin}
                  \item{\code{Rel}}{semantic similarty of Schlicker \emph{et al.}}
                  \item{\code{JC}}{semantic similarty of Jiang and Conrath}
                  \item{\code{Nunivers}}{semantic similarty of Mazandu and Mulder}
                  \item{\code{Psec}}{semantic similarty of Pirro and Seco}
                  \item{\code{Faith}}{semantic similarty of Pirro and Euzenat}
                  \item{\code{all}}{all sementic similarities}
                }
             By default computes the semantic similarity of Resnik, that is, \code{Res}}
}
\details{This function computes either all the semantic similarities implemented in the packages or the one indicated.}
\value{The resulting object is an \code{data.frame} where for each pair of different terms (rows) are shown the \code{numeric} values of the semantic similarities computed (columns)}
%\references{}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{ancestors}}, \code{\link{simRes}}, \code{\link{simLin}}, \code{\link{simRel}}, \code{\link{simJC}}, \code{\link{simNunivers}}, \code{\link{simPsec}}, \code{\link{simFaith}}, \code{\link{summaryMICA}}, \code{\link{sims.eb}}, \code{\link{pseudoDists}},\code{\link{resnikSummary}}}
\examples{
data(joslyn)

inv.IminusG <- inverseIminusG(joslyn.OOC)
A.mat <- getA(inv.IminusG)
at <- ancestors(A.mat)
resnik.sum <- resnikSummary(x = joslyn.OOC)
ic <- resnik.sum[, "ic"]

sims.nb(at, ic, method = "Res")     # Renik
sims.nb(at, ic, method = "Lin")     # Lin
sims.nb(at, ic, method = "Rel")     # Schlicker et al.
sims.nb(at, ic, method = "JC")      # Jiang and Conrath
sims.nb(at, ic, method = "Nuviers") # Mazandu and Mulder
sims.nb(at, ic, method = "Psec")    # Pirro & Seco
sims.nb(at, ic, method = "Faith")   # Pirro & Euzenat

sims.nb(at, ic, method = "all")     # All semantic similarities
}
\keyword{ sims.nb }
\keyword{ sims.eb }
\keyword{ pseudoDists }
\keyword{ ancestors }
\keyword{ summaryMICA }
\keyword{ simRes }
\keyword{ simLin }
\keyword{ simtJC }
\keyword{ simNunivers }
\keyword{ simPsec }
\keyword{ simFaith }
\keyword{ resnikSummary }