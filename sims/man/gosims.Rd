\name{gosims}
\alias{gosims}
\title{Wrapper function that calls different approaches and methods for computing semantic similarities between GO Identifiers given a list of Entrez Gene IDs}
\description{Given a list of Entrez Gene Identifiers, the name of an \code{R} organism package, and an ontology domain, the function builds an Object-Ontology Complex particularly adapted for dealing with GO Identifiers. But also, if two list of Entrez Gene Identifiers are provided, then the function builds a list with two elements, each one associated with one of the list of Entrez Gene Identifiers and consisting of an object of class \code{\link{OOC}}}
\usage{gosims(eg, ontology = "BP", pckg = "org.Hs.eg.db", type = "nb", method = "Res")
}
\arguments{
	\item{eg}{\code{character} vector with the Entrez Gene IDs to be interrogated.}
	\item{ontology}{\code{character} with the ontology domain of the GO (see \emph{Details}). By default is \code{BP}}
	\item{pckg}{\code{character} with the name of the \code{R} organism package. By default is \code{org.Hs.eg.db}}
	\item{type}{\code{character} indicating the approach for computing semantic similarities. By default \code{nb}}
	\item{method}{\code{character} indicating the method used for computing semantic similarities. By default \code{Res}}
}
\details{
Argument \code{ontology} has three possible domains
\describe{
	\item{\code{BP}}{Biological Processes}
	\item{\code{CC}}{Cellular Components}
	\item{\code{MF}}{Molecular Functions}
}

Argument \code{type} has three possibilities
\describe{
	\item{\code{nb}}{node-based approach}
	\item{\code{eb}}{edge-based approach}
	\item{\code{eb.pd}}{edge-based approach based on pseudo-distances}
}
Depending on the approach selected in argument \code{type}, different \code{methods} for computing semantic similarities are available.
\itemize{
	\item{Possible methods of \emph{Node-based approach} (i.e. \code{type = 'nb'}) are
		\describe{
        		\item{\code{Res}}{semantic similarty of Resnik}
        		\item{\code{Lin}}{semantic similarty of Lin}
        		\item{\code{Rel}}{semantic similarty of Schlicker \emph{et al.}}
        		\item{\code{JC}}{semantic similarty of Jiang and Conrath}
        		\item{\code{Nunivers}}{semantic similarty of Mazandu and Mulder}
        		\item{\code{Psec}}{semantic similarty of Pirro and Seco}
        		\item{\code{Faith}}{semantic similarty of Pirro and Euzenat}
        		\item{\code{all}}{all sementic similarities}
        	}}
	\item{Possible methods of \emph{Edge-based approach} (i.e. \code{type = 'eb'}) are
		\describe{
        		\item{\code{Rada}}{semantic similarty of Rada}
            \item{\code{Res.eb}}{semantic similarty of Resnik based on the shortest path (see function \code{\link{distRada}})}
            \item{\code{all}}{all semantic similarities of edge-based approaches implemented in the package}
       	}}
    	\item{Possible pseudo-distance methods of \emph{Edge-based approach} (i.e. \code{type = 'eb.pd'}) are
		\describe{
         	\item{\code{hm}}{pseudo-distance of the minimum chain length}
            \item{\code{hx}}{pseudo-distance of the maximum chain length}
            \item{\code{hax}}{pseudo-distance of the average of extreme chain lengths}
            \item{\code{hap}}{pseudo-distance of the average of all chain lengths}
            \item{\code{all}}{all pseudo-distances}
		}}
}
}
%\value{}
%\references{}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{gosimsAvsB}}, \code{\link{simsBetweenGOIDs}}, \code{\link{summarySims}}, \code{\link{goOOC}}, \code{\link{sims.nb}}, \code{\link{sims.eb}}, \code{\link{pseudoDists}}}
\examples{
data(prostateIds) 

## Entrez Genes from two different studies
## of prostate cancer
eg.we <- welsh01EntrezIDs[1:10]   # Welsh study
eg.sg <- singh01EntrezIDs[1:10]   # Singh study

## All semantic similarities from node-based approach

nb.all <- gosims(eg = eg.we, ontology = "MF", 
                 pckg = "org.Hs.eg.db", type = "nb", 
                 method = "all")
head(nb.all, 20)
tail(nb.all, 20)

## All semantic similarities from edge-based approach

eb.all <- gosims(eg = eg.we, ontology = "MF", 
                 pckg = "org.Hs.eg.db", type = "eb", 
                 method = "all")
head(eb.all, 20)
tail(eb.all, 20)

## All pseudo-distances from node-based approach

pd.all <- gosims(eg = eg.we, ontology = "MF", 
                 pckg = "org.Hs.eg.db", 
                 type = "eb.pd", method = "all")
head(pd.all, 20)
tail(pd.all, 20)
}
\keyword{ gosims }
\keyword{ gosimsAvsB }
\keyword{ simsBetweenGOIDs }
\keyword{ summarySims }
\keyword{ goOOC }
\keyword{ sims.nb }
\keyword{ sims.eb }
\keyword{ pseudoDists }
