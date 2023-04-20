\name{toOOC}
\alias{toOOC}
\title{Builds an Object-Ontology Complex (OOC)}
\description{Builds a object of class \code{OOC} that is used as a container of an Object-Ontology Object (OOC)}
\usage{toOOC(T, G, O, M)}
\arguments{
    \item{T}{\code{character} vector with the names of the terms in the ontology}
    \item{G}{\code{matrix} encoding the refinement matrix associated with the DAG of the ontology. Elements of this 
                    matrix are 1, when a term on the rows refines a term on the columns, and 0 otherwise}
    \item{O}{\code{character} vector with the object identifiers (i.e. features, genes,...) mapping the terms of the
                    ontology}
    \item{M}{\code{matrix} encoding the mapping from the set of object identifiers (listed in the argument \code{O} to 
                    the terms listed in the argument \code{T}}
}
%\details{}
\value{The resulting object is an \code{OOC}}
%\references{}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{OOC}}}
\examples{
## Terms in the ontology

terms <- c("R", "B", "C", "K", "F", "G", "I", "E", "J", "H", "A", "D")

## Object identifiers

object.ids <- letters[1:10]

## Matrix of refinements

G <- matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,
              1,0,0,0,0,0,0,0,0,0,0,0,
              1,0,0,0,0,0,0,0,0,0,0,0,
              1,0,0,0,0,0,0,0,0,0,0,0,
              0,1,0,0,0,0,0,0,0,0,0,0,
              0,1,0,0,0,0,0,0,0,0,0,0,
              0,1,1,0,0,0,0,0,0,0,0,0,
              0,0,1,1,0,0,1,0,0,0,0,0,
              0,0,1,1,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,1,0,0,0,0,0,
              0,0,0,0,1,1,0,0,0,1,0,0,
              0,0,0,0,0,0,0,1,1,0,0,0),
             nrow = 12, ncol = 12, byrow = TRUE)
colnames(G) <-rownames(G) <- terms

## Mapping matrix

M <- matrix(c(0,0,0,0,0,0,0,0,0,0,1,0,
              0,0,0,0,1,0,0,1,0,0,1,0,
              0,0,0,0,0,0,0,0,0,0,1,0,
              0,0,0,0,1,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,1,0,0,
              0,0,0,0,0,0,1,0,0,0,0,0,
              0,0,0,0,0,0,0,0,1,0,0,0,
              0,0,0,0,0,0,0,0,1,0,0,0,
              0,0,0,0,0,0,0,0,1,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,1),
             nrow = 10, ncol = 12, byrow = TRUE)
colnames(M) <- terms
rownames(M) <- object.ids

## Object-Ontology Complex

joslyn.OOC <- toOOC(T = terms, G = G, O = object.ids, M = M)
}
\keyword{ toOOC }
\keyword{ OOC }
