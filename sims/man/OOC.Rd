\name{OOC}
\Rdversion{1.1}
\docType{class}
\alias{OOC}
\title{General container for an Object-Ontology Complex (OOC)}
\description{Class \code{"OOC"} is a general container for Object-Ontology Complexes (OOC)}
\section{Objects from the Class}{
Objects can be created by calls of either the form \code{new("OOC", ...)} or the function \code{toOOC(T, G, O, M)}}
\section{Slots}{
  \describe{
    \item{\code{T}}{\code{character} vector with the names of the terms in the ontology}
    \item{\code{G}}{\code{matrix} encoding the refinement matrix associated with the DAG of the ontology. Elements of this 
                    matrix are 1, when a term on the rows refines a term on the columns, and 0 otherwise}
    \item{\code{O}}{\code{character} vector with the object identifiers (i.e. features, genes,...) mapping the terms of the
                    ontology}
    \item{\code{M}}{\code{matrix} encoding the mapping from the set of object identifiers (listed in the slot \code{O} to the 
                     terms listed in the slot \code{T}}}
}

\section{Methods}{No methods defined with class \code{OOC} in the signature}
%\references{}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{toOOC}}}
\keyword{ OOC }
\keyword{ toOOC }