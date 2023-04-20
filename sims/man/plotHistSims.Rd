\name{plotHistSims}
\alias{plotHistSims}
\title{Histogram of two semantic similarity profiles}
\description{This function plots an histogram of two semantic similarity profiles in the same plot}
\usage{plotHistSims(x, freq = TRUE, main = "Histogram of Semantic Similarities", xlab = "Semantic Similarity")}
\arguments{
	\item{x}{\code{matrix} with measures of semantic similarity measures provided by the function \code{\link{gosimsAvsB}}}
	\item{freq}{\code{logical} value. If \code{TRUE}, histogram graphic is a representation of frequencies, the ``counts'' 
	            component of the result. If \code{FALSE}, probability densities, component ``density'', are plotted (so that 
	            the histogram has a total area of one). By default is \code{TRUE}}
	\item{main}{\code{character} to title the graphic}
	\item{xlab}{\code{character} to label x-axis}
}
\details{The plot shows two ``curves'' in the same figure overlapped. Area colored in \code{red} is associated with the first column of the \code{matrix} provided in argument \code{x}, and area colored in \code{blue} is associated with the second column of the \code{matrix}. Intersected area is colored in \emph{violet}.}
\value{Resulting object is of class \code{histogram}}
%\references{}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{gosimsAvsB}}, \code{\link{summarySimsAvsB}}, \code{\link{gosimsProfiles}}}
\examples{
data(prostateIds)                # Data set from the package goProfiles

pckg = "org.Hs.eg.db"            # Organism package of humans

eg.we <- welsh01EntrezIDs[1:10]   # Entrez Gene IDs signature Welsh 01
eg.sg <- singh01EntrezIDs[1:10]   # Entrez Gene IDs signature Singh 01

WEvsSG.nb <- gosimsAvsB(eg1 = eg.we, eg2 = eg.sg, 
                        ontology = "MF", pckg = pckg, 
                        type = "nb", method = "Res")

plotHistSims(x = WEvsSG.nb, freq = TRUE, 
             main = "Histogram of Semantic Similarities",
             xlab = "Semantic Similarity")


}
\keyword{ plotHistSims }
\keyword{ plotGODAG }
\keyword{ gosimsProfiles }
\keyword{ gosimsAvsB }
\keyword{ summarySimsAvsB }

