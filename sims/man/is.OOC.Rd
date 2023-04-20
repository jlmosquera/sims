\name{is.OOC}
\alias{is.OOC}
\title{Tests if its argument is a (strict) \code{OOC} object}
\description{Tests if its argument is a (strict) \code{OOC} object}
\usage{is.OOC(x)}
\arguments{
  \item{x}{an R object}
}
\details{\code{is.OOC} returns \code{TRUE} if argument \code{x} is an \code{OOC} object, and \code{FALSE} otherwise}
%\value{}
%\references{}
\author{Jose Luis Mosquera}
%\note{}
\seealso{\code{\link{toOOC}}}
\examples{
data(joslyn)

is.OOC(mat.m)            # FALSE
is.OOC(joslyn.OOC)   # TRUE
}
\keyword{ is.OOC }
\keyword{ toOOC }
