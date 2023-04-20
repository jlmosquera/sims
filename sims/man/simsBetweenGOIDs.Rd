\name{simsBetweenGOIDs}
\alias{simsBetweenGOIDs}
\title{Wrapper function that calls different approaches and methods for computing semantic similarities between GO ID ancestors of a list of GO ID's}
\description{Given a list of GO ID's, a GO ontology domain, a type of approach for computing semantic similarities and a specific method, the function calculates the semantic similarities between all the pairs of GO ID ancestors of the induced graph.}
\usage{simsBetweenGOIDs(goids, ontology = "BP", type = "nb", method = "Res")
}
\arguments{
	\item{goids}{\code{character} vector with the GO IDs.}
	\item{ontology}{\code{character} with the ontology domain (see \emph{Details}). By default is \code{BP}}
	\item{type}{\code{character} indicating the approach for computing semantic similarities. By default \code{nb}}
	\item{method}{\code{character} indicating the method used for computing semantic similarities. By default \code{Res}}
}
\details{
There are three possibilities, of course, in argument \code{ontology}
		\describe{
			\item{\code{BP}}{Biological Processes}
			\item{\code{CC}}{Cellular Components}
			\item{\code{MF}}{Molecular Functions}
		}

There are three possibilities in argument \code{type}
		\describe{
			\item{\code{nb}}{node-based approach}
			\item{\code{eb}}{edge-based approach}
			\item{\code{eb.pd}}{edge-based approach based on psedudo-distances}
		}

According to the approach selected in argument \code{type}, there are different \code{methods} for computing semantic similarities.
\itemize{
	\item{Node-based approach (i.e. \code{type = 'nb'})
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
	\item{Edge-based approach (i.e. \code{type = 'eb'})
		\describe{
        		\item{\code{Rada}}{semantic similarty of Rada}
            \item{\code{Res.eb}}{semantic similarty of Resnik based on the shortest path (see function \code{\link{distRada}}}
            \item{\code{all}}{all semantic similarities of edge-based approaches implemented in the package}
       	}}
    	\item{Edge-based approach pseudo-distance (i.e. \code{type = 'eb.pd'})
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
\seealso{\code{\link{gosims}}, \code{\link{gosimsAvsB}}}
\examples{
goids <- c("GO:0004022","GO:0005515")
a <- simsBetweenGOIDs(goids, ontology = "MF", type = "nb", method = "Res")
}
\keyword{ simsBetweenGOIDs }
\keyword{ gosims }
\keyword{ gosimsAvsB }
