\name{goOOC}
\alias{goOOC}
\title{Builds an Object-Ontology Complex (\code{OOC}) whose slots are associated with GO Identifiers}
\description{Given a list of Entrez Gene Identifiers, the name of an \code{R} organism package, and a domain of the GO, the function builds an Object-Ontology Complex where the slots are associated with GO ID's and Entrez Gene ID's. But also, if two list of Entrez Gene Identifiers are provided, then the function builds a list with two elements, each one associated with one of the list of Entrez Gene Identifiers and consisting of an object of class \code{\link{OOC}}}
\usage{goOOC(eg1, eg2 = NULL, pckg = "org.Hs.eg.db", ontology = "BP")}
\arguments{
	\item{eg1}{\code{character} vector with the (first list of) Entrez Gene IDs to be interrogated.}
	\item{eg2}{\code{character} vector with the second list of Entrez Gene IDs to be interrogated. By default is \code{NULL}.}
	\item{pckg}{\code{character} with the name of the \code{R} organism package. By default is \code{org.Hs.eg.db}}
	\item{ontology}{\code{character} the ontology from the GO is selected (see \emph{Details}. By default is \code{BP}}
}
\details{This function is similar to the function \code{\link{toOOC}}, but it is particularly adapted to deal with information associated with the Gene Ontology.

If argument \code{eg2} is \code{NULL}, then yields an \code{\link{OOC}} object. Otherwise, the resulting object is a \code{list} with two \code{\link{OOC}} objects. First element is associated with the first list of Entrez Gene Identifiers (provided to the function in argument \code{eg1}), and second element is associated with the second list of Entrez Gene Identifiers (provided to the function in argument \code{eg2}).

Note that, in case of providing two lists of Entrez Gene Identifiers, both objects \code{\link{OOC}} from the resulting list have the same refinement matrix (i.e. slot \code{G}). This is made in order to compare both lists of Entrez Gene IDs in terms of semantic similarities profiles later on. That is, sometimes the comparison between two list of genes is required. This tasks may be performed by comparing two semantic similarity profiles (i.e. two lists of semantic similarities between each pair of terms and performed with the same measure). Thus, each list of genes yields an induced subgraph from the GO. But, note that both subgraphs might be different. In consequence, no comparison might be computed because the number of pairs of GO IDs wolud be different. This drawback, can be solved by considering a common subgraph induced by both lists of genes.

Argument \code{ontology} has three possibilities
\describe{
	\item{\code{BP}}{Biological Processes}
	\item{\code{CC}}{Cellular Components}
	\item{\code{MF}}{Molecular Functions}
}
}
\value{Resulting object}
%\references{}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{toOOC}}}
\examples{
data(prostateIds)

## Entrez Genes from two different studies
## of prostate cancer
eg.we <- welsh01EntrezIDs[1:10]  # Weslh study
eg.sg <- singh01EntrezIDs[1:10]  # Singh study

## OOC associated with Entrez Gene IDs of Welsh

ooc1 <- goOOC(eg1 = eg.we, eg2 = NULL, 
              pckg = "org.Hs.eg.db", ontology = "MF")
class(ooc1)
str(ooc1)
dim(ooc1@G)
dim(ooc1@M)

## OOC associated with Entrez Gene IDs of Singh
        
ooc2 <- goOOC(eg1 = eg.sg, eg2 = NULL, 
              pckg = "org.Hs.eg.db", ontology = "MF")
class(ooc2)
str(ooc2)
dim(ooc2@G)
dim(ooc2@M)

## List of OOC associated with Entrez Gene IDs of Welsh 
## and Singh
        
ooc <- goOOC(eg1 = eg.we, eg2 = eg.sg, 
             pckg = "org.Hs.eg.db", ontology = "MF")
class(ooc)
str(ooc)
dim(ooc[[1]]@G)
dim(ooc[[2]]@G)
table(ooc[[1]]@G==ooc[[2]]@G)    
dim(ooc[[1]]@M)
dim(ooc[[2]]@M)
}
\keyword{ goOOC }
\keyword{ toOOC }
