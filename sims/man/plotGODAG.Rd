\name{plotGODAG}
\alias{plotGODAG}
\title{Plots a subgraph from the GO associated with one or two lists of Entrez Gene Identifiers}
\description{This function plots a subgraph from the ontology selected from the GO (see \emph{Details}) associated with one or two lists of Entrez Gene IDs}
\usage{plotGODAG(eg1, eg2 = NULL, pckg = "org.Hs.eg.db", ontology = "MF", verbose = FALSE)}
\arguments{
	\item{eg1}{\code{character} vector with the (first list of) Entrez Gene IDs}
	\item{eg2}{\code{character} vector with the second list of Entrez Gene IDs. By default is \code{NULL}}
	\item{pckg}{\code{character} with the name of the \code{R} organism package. By default is \code{org.Hs.eg.db}}
	\item{ontology}{\code{character} the ontology from the GO is selected (see \emph{Details}. By default is \code{BP}}
	\item{verbose}{\code{logical} . By default \code{FALSE}}
}
\details{
There are three possibilities, of course, in argument \code{ontology}
		\itemize{
			\item{BP}{Biological Processes}
			\item{CC}{Cellular Components}
			\item{MF}{Molecular Functions}
		}

The subgraph shows two types of shapes for each node. Circles are nodes not mapped directly and rectangles are nodes mapped directly.

The color of nodes indicate the type of relation with the Entrez Gene IDs. That is, when argument \code{eg2} is \code{NULL}, there are two possibilities: nodes mapped directly are shown in  \code{red} color and their ancestors are shown in \code{yellow} color. But, if argument \code{eg2} is not \code{NULL}, then there are six different colors. Nodes mapped directly from the first list of Entrez Gene IDs are shown in \code{red} color and their ancestors are shown in \code{yellow} color. Nodes mapped directly from the second list of Entrez Gene IDs are shown in \code{lightblue} color and their ancestors are shown in \code{blue} color. Nodes mapped directly from both lists of Entrez Gene IDs are shown in \code{magenta} color and their ancestors are shown in \code{violet} color.}
%\value{}
%\references{}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{gosimsAvsB}}, \code{\link{summarySimsAvsB}}, \code{\link{gosimsProfiles}}}
\examples{
data(prostateIds)                # Data set from the package goProfiles

pckg = "org.Hs.eg.db"            # Organism package of humans

eg.we <- welsh01EntrezIDs[1:10]   # Entrez Gene IDs signature Welsh 01
eg.sg <- singh01EntrezIDs[1:10]   # Entrez Gene IDs signature Singh 01

plotGODAG(eg1 = eg.we, eg2 = NULL, pckg = pckg, ontology = "MF")
plotGODAG(eg1 = eg.sg, eg2 = NULL, pckg = pckg, ontology = "MF")
plotGODAG(eg1 = eg.we, eg2 = eg.sg, pckg = pckg, ontology = "MF")
}
\keyword{ plotGODAG }
\keyword{ gosimsProfiles }
\keyword{ gosimsAvsB }
\keyword{ summarySimsAvsB }

