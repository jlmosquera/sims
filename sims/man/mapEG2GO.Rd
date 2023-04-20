\name{mapEG2GO}
\alias{mapEG2GO}
\title{Mapping Entrez Gene IDs to Gene Ontology IDs}
\description{This function maps the Entrez Gene IDs on to the GO IDs provided by an \code{R} organism package (e.g. \code{org.Hs.eg.db}}
\usage{mapEG2GO(eg = NULL, pckg = "org.Hs.eg.db")}
\arguments{
    \item{eg}{\code{character} vector with the Entrez Gene IDs. By default is \code{NULL}, which means that all the eg from the \code{pckg} will be mapped}
    \item{pckg}{\code{character} with the name of the package of the organism from where the GO IDs are extracted. By default is \code{org.Hs.eg.db}, that is, Homo sapiens specie}
}
\value{The resulting object is a \code{list} with three elements, corresponding to an ontology of the the GO (i.e \code{BP}, \code{CC} and \code{MF}). Each of these elements is a \code{data.frame} with two columns. First column consists of the list of Entrez Gene IDs mapping the GO IDs placed in the second column. Note that in both columns, the concepts can be repeated. That is, an Entrez Gene ID can map different GO IDs, and a GO ID can be mapped by different Entrez Gene IDs}
\author{Jose Luis Mosquera}
%\note{}
%\seealso{\code{\link{}}, \code{\link{}}}
\examples{
data(prostateIds)               # Data set from the package goProfiles

pckg <- "org.Hs.eg.db"          # Names of the organism package
eg <- welsh01EntrezIDs[1:10]    # Entrez Gene IDs

eg2go <- mapEG2GO(eg, pckg)
lapply(eg2go, head)
}
\keyword{ mapEG2GO }

