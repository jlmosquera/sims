\name{toMat}
\alias{toMat}
\title{Builds a \code{matrix} of zero and one elements such that zero indicates there is no a relation between row and column, and one there is a relation.}
\description{Givena two \code{character} vectors and a 2-columns \code{data.frame}, that relates an element of the first \code{character} vector with an element of the second \code{character} vector, builds a \code{matrix} whose \code{rownames} are the elements of the first \code{character} vector, \code{colnames} are the elements of the second \code{character} vector, and elements are 0 and 1. The values indicate that there is no a relation between a row and a column, when the value is zero, and there is a relation between a row and a column, when the value is one}
\usage{toMat(df, rnames, cnames)}
\arguments{
    \item{df}{\code{data.frame} with mapping from the Entrez Gene IDs to the GO ID terms}
    \item{rnames}{\code{character} indicating which ontologies is selected (see \emph{Details})}
    \item{cnames}{\code{character} indicating which ontologies is selected (see \emph{Details})}
}
%\details{}
\value{The resulting object is a \code{matrix}}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{toPairs}}, \code{\link{mappingMatrix}}, \code{\link{refinementMatrix}}}
\examples{
vocabulary <- c("R", "B", "C", "K", "F", "G", "I", "E", "J", "H", "A", "D")
object.ids <- letters[1:10]

origin <- c("B", "C", "K", "F", "G", "I", "I", "E", "J", "E", "J", "A", "A",
            "E", "H", "D", "D", "A")
terminus <- c("R", "R", "R", "B", "B", "B", "C", "C", "C", "K", "K", "F", "G",
              "I", "I", "E", "J", "H")
links <- data.frame(origin, terminus)

mat.g <- toMat(df = links, rnames = vocabulary, cnames = vocabulary)
}
\keyword{ toMat }
\keyword{ toPairs }
\keyword{ refinementMatrix }
\keyword{ mappingMatrix }