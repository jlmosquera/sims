\name{sims}
\alias{sims}
\docType{package}
\title{Computation of semantic similarity measures between terms of an ontology}
\description{\code{sims} package is devoted to compute semantic similarity measures between terms of an ontology, and some functions are particularly addressed to manage the Gene Ontology terms.}

\details{This package is designed to compute semantic similarities between terms of any ontology. Fourteen measures from different approaches are implemented. Specifically, from node-based approach there are implemented seven semantic similarity measures proposed by Resnik, Lin, Schlicker et al., Jiang and Conrath, Mazandu and Mulder, Pirro and Seco, and Pirro and Euzenat. With regard to edge-based approaches there are implemented: two semantic similarity measures proposed by Resnik, and Rada et al, one distance measure proposed by Rada and four pseudo-distances proposed by Joslyn et al. The package can manage any ontology, but it is particularly focused on the Gene Ontology. In this regard, there are some functions that allow building the refinements matrix (i.e. the accessibility matrix in terms of graph theory), the mapping matrix (i.e the matrix that maps from Entrez Gene IDs to GO IDs) and plot the DAG structure in order to compare two different list of GO terms. \code{sims} package can manage Entrez Gene IDs and GO IDs from any organism \code{R} package. 

\tabular{ll}{
Package: \tab sims\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 31-08-2014\cr
License: \tab GPL-2 \cr
Depends: \tab AnnotationDbi, expm, goProfiles, GOstats, igraph, methods, plotrix, Rgraphviz, vegan\cr
Imports: \tab Matrix, plyr, knitr\cr
Suggests: \tab org.Hs.eg.db\cr
}
The most important functions are \code{sims.nb}, \code{sims.eb}, \code{pseudoDists}, \code{gosims} and \code{gosimsAvsB}
}
\author{
Jose Luis Mosquera and Alex Sanchez

Maintainer: Jose Luis Mosquera <jlmosquera@gmail.com>
}
%\references{}
\seealso{\code{\link{sims.nb}}, \code{\link{sims.eb}}, \code{\link{pseudoDists}}, \code{\link{gosims}}, \code{\link{gosimsAvsB}}}
\keyword{ package }
\keyword{ sims }