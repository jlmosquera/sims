\name{refinementMatrix}
\alias{refinementMatrix}
\title{Builds the refinement matrix associated with the DAG structure of Gene Ontology}
\description{Based on the \code{data.frame} of the maps the Entrez Gene IDs on to the GO IDs, this function builds the refinement matrix (i.e. in terms of the graph theory the accessibility matrix) associated with the GO DAG structure}
\usage{refinementMatrix(df, ontology = "BP")}
\arguments{
    \item{df}{\code{data.frame} with mapping from the Entrez Gene IDs to the GO ID terms}
    \item{ontology}{\code{character} indicating which ontologies is selected (see \emph{Details}). By default is \code{BP} (i.e. Biological Processes}
}
\details{Given the \code{data.frame} with mapping from the Entrez Gene IDs to the GO IDs, \code{refinementMatrix} looks for all the ancestors and builds the \code{matrix} with the refinements associated with the directed subgraph extracted from the GO DAG.

\code{ontology} argument requires to indicate which of the ontologies in GO is selected for building the refinement matrix. Obviously, there are three possibilities
\itemize{
	\item{BP}{Biological Processes}
	\item{CC}{Cellular Components}
	\item{BP}{Molecular Functions}		
}
}
\value{The resulting object is a binary \code{matrix} whose rows and columns are all the GO IDs that make up the subgraph of the GO DAG. Each element of this \code{matrix} is either 1 when there exists a refinement from the GO ID on the row to the GO ID on the column, or 0 otherwise.}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{mapEG2GO}}, \code{\link{mappingMatrix}}}
\examples{
data(prostateIds)               # Data set from the package goProfiles

pckg <- "org.Hs.eg.db"          # Names of the organism package

eg <- welsh01EntrezIDs[1:10]    # Entrez Gene IDs
eg2go <- mapEG2GO(eg, pckg)
eg2go.mf <- eg2go$MF

G <- refinementMatrix(df = eg2go.mf, ontology = "MF")
print(G)
}
\keyword{ refinementMatrix }
\keyword{ mapEG2GO }
\keyword{ mappingMatrix }