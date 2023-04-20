\name{mappingMatrix}
\alias{mappingMatrix}
\title{Mapping matrix from the Entrez Gene IDs to the GO IDs associated with the directed subgraph extracted from GO DAG structure}
\description{Based on the \code{data.frame} of the maps the Entrez Gene IDs on to the GO IDs and the refinement matrix, this function builds the mapping matrix associated with the GO DAG structure}
\usage{mappingMatrix(G, df)}
\arguments{
	\item{G}{binary \code{matrix} indicating the refinements of each GO ID}
    \item{df}{\code{data.frame} with mapping from the Entrez Gene IDs to the GO ID terms}
}
\details{Given the \code{data.frame} with mapping from the Entrez Gene IDs to the GO IDs, \code{mappingMatrix} and the refinement matrix, builds a binary \code{matrix} where on rows are the Entrez Gene IDs and on the columns are all the GO IDs associated with the directed subgraph extracted from the GO DAG.}
\value{The resulting object is a binary \code{matrix}. Each element of this \code{matrix} is either 1 when an Entrez Gene ID on a row is mapping a GO ID on a column, or 0 otherwise.}
\author{Jose Luis Mosquera}
%\note{
%%  ~~further notes~~
%}
\seealso{\code{\link{mapEG2GO}}, \code{\link{refinementMatrix}}}
\examples{
data(prostateIds)                     # Data set from the package goProfiles

pckg <- "org.Hs.eg.db"                # Names of the organism package
eg <- welsh01EntrezIDs[1:10]          # Entrez Gene IDs

eg2go <- mapEG2GO(eg, pckg)
eg2go.mf <- eg2go$MF                  # Get mollecular functions

G <- refinementMatrix(df = eg2go.mf, ontology = "MF")
M <- mappingMatrix(G, df = eg2go.mf)
print(M)
}
\keyword{ mappingMatrix }
\keyword{ mapEG2GO }
\keyword{ refinementMatrix }
