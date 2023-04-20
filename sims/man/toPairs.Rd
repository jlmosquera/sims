\name{toPairs}
\alias{toPairs}
\title{Builds a 2-columns \code{data.frame} relating elements of a matrix with value one.}
\description{Given a matrix of zero's and one's, builds a 2-columns \code{data.frame} where on the first columns are the names of the rows and on the second column the names of columns for those elements of the matrix with value one.}
\usage{toPairs(mat)}
\arguments{
    \item{mat}{\code{matrix} with 0's and 1's}
}
%\details{}
\value{The resulting object is a 2-columns \code{data.frame}}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{toMat}}, \code{\link{mappingMatrix}}, \code{\link{refinementMatrix}}}
\examples{
a <- sample(0:1, 20, replace = TRUE)
A <- matrix(a, nrow = 5)
colnames(A) <- LETTERS[1:4]
rownames(A) <- letters[1:5]

toPairs(A)
}
\keyword{ toPairs }
\keyword{ toMat }
\keyword{ refinementMatrix }
\keyword{ mappingMatrix }