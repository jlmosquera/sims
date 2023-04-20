\name{gosimsAvsB}
\alias{gosimsAvsB}
\title{Wrapper function for computing semantic similarities between GO Identifiers for two lists of Entrez Gene IDs}
\description{Given two lists of Entrez Gene Identifiers computes semantic similarities between GO Identifiers}
\usage{gosimsAvsB(eg1, eg2, ontology = "BP", pckg = "org.Hs.eg.db", type = "nb", method = "Res")
}
\arguments{
	\item{eg1}{\code{character} vector with the first list of Entrez Gene IDs to be interrogated.}
	\item{eg2}{\code{character} vector with the second list of Entrez Gene IDs to be interrogated.}
	\item{ontology}{\code{character} the ontology from the GO is selected (see \emph{Details}). By default is \code{BP}}
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
\value{Resulting object is a \code{matrix} semantic similarity values. On the rows there are the pairs of GO ID's compared and on the columns the method(s) used for the computation of the measures. For each measure there are two columns (one per list of genes)}
%\references{}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{summarySims}}, \code{\link{gosims}}, \code{\link{goOOC}}, \code{\link{sims.nb}}, \code{\link{sims.eb}}, \code{\link{pseudoDists}}}
\examples{
data(prostateIds)

## R organism package
pckg <- "org.Hs.eg.db"

## Entrez Genes from two different studies
## of prostate cancer
eg.we <- welsh01EntrezIDs[1:10]   # Welsh study
eg.sg <- singh01EntrezIDs[1:10]   # Singh study

## Semantic similarity of Resnik (node-based approach)

WEvsSG.nb <- gosimsAvsB(eg1 = eg.we, eg2 = eg.sg, 
                        ontology = "MF", pckg = pckg, 
                        type = "nb", method = "Res")
head(WEvsSG.nb, 20)
tail(WEvsSG.nb, 20)

## Semantic similarity of Resnik (edge-based approach)

WEvsSG.eb <- gosimsAvsB(eg1 = eg.we, eg2 = eg.sg, 
                        ontology = "MF", pckg = pckg, 
                        type = "eb", method = "Res.eb")
head(WEvsSG.eb, 20)
tail(WEvsSG.eb, 20)

## Pseudo-distance of the minimum chains length 
## (edge-based approach)

WEvsSG.pd <- gosimsAvsB(eg1 = eg.we, eg2 = eg.sg, 
                        ontology = "MF", pckg = pckg, 
                        type = "eb.pd", method = "hm")
head(WEvsSG.pd, 20)
tail(WEvsSG.pd, 20)
}
\keyword{ gosimsAvsB }
\keyword{ summarySims }
\keyword{ gosims }
\keyword{ goOOC }
\keyword{ sims.nb }
\keyword{ sims.eb }
\keyword{ pseudoDists }
