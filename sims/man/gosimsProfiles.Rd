\name{gosimsProfiles}
\alias{gosimsProfiles}
\title{Plots a vertical bar diagram whose bars are associated with the semantic similarities between each pair of terms, and such that bars on the left side of the plot are the corresponding to the first group of objects and on the bars on the right side are the bars corresponding to the second group of objects}
\description{Given a two-columns \code{matrix} (or \code{data.frame}), this function yields a vertical bar diagram whose bars are associated with the semantic similarities between each pair of terms, and such that bars on the left side of the plot are the corresponding to the first group of objects and on the bars on the right side are the bars corresponding to the second group of objects}
\usage{gosimsProfiles(x, col = c("tomato", "blue"), cex = 0.4, 
                      top.labels = NULL, main = NULL, xlab = "Semantinc Similarity")}
\arguments{
	\item{x}{\code{matrix} (or \code{data.frame}) where for each pair terms (rows) contains the values of the semantic similarity measure estimated for each group of objects (columns)}
	\item{col}{\code{character} vector with the two elements containing the colors associated with each group of objects. By default are \code{tomato} for the left bars on the left side, and \code{blue} for the bars on the right side}
	\item{cex}{\code{numeric} value with the expansion for the category labels. By default is 0.4}
	\item{top.labels}{\code{character} vector indicating the two categories represented on the left and right sides of the plot and a heading for the labels in the center. By default is \code{NULL}}
	\item{main}{\code{character} with title for the plot. By default \code{NULL}}
	\item{xlab}{\code{character} with the label for the units of the plot (i.e. labels on the x-axis).}
}
%\details{}
%\value{}
%\references{}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{gosimsAvsB}}, \code{\link{summarySimsAvsB}}, \code{\link{plotHistSims}}, \code{\link{plotGODAG}}}
\examples{
data(prostateIds)

## Entrez Genes from two different studies
## of prostate cancer
eg.we <- welsh01EntrezIDs[1:10]   # Welsh study
eg.sg <- singh01EntrezIDs[1:10]   # Singh study

WEvsSG.nb <- gosimsAvsB(eg1 = eg.we, eg2 = eg.sg, 
                        ontology = "MF", pckg = "org.Hs.eg.db",
                        type = "nb", method = "Res")

## Plot of the whole profiles of semantic similarities

gosimsProfiles(x = WEvsSG.nb, col = c("tomato", "blue"), 
               cex = 0.4, top.labels = c("WE", "SG"),
               main = "Welsh 01 vs Singh 01", xlab = "Resnik")

## Plot a subset of the whole profiles of semantic similarities

gosimsProfiles(x = WEvsSG.nb[1:200, ], col = c("tomato", "blue"), 
               cex = 0.4, top.labels = c("WE", "SG"),
               main = "Welsh 01 vs Singh 01", xlab = "Resnik")
}
\keyword{ gosimsProfiles }
\keyword{ gosimsAvsB }
\keyword{ summarySimsAvsB }
\keyword{ plotHistSims }
\keyword{ plotGODAG }
