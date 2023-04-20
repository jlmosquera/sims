\name{Nt}
\alias{Nt}
\title{Number of times that each term or any of its specializations references to an ancestor}
\description{Calulates the \code{matrix} of the number of times that each term or any of its specializations references to an ancestor}
\usage{Nt(x)}
\arguments{
  \item{x}{an object of class \code{OOC}}
}
\details{The resulting matrix is based on the Neumann series of a matrix. The function calculates the \code{matrix} based on the formula 
\deqn{N_{t}=M(I-\Gamma)^{-1}}
}
%\value{}
%\references{}
\author{Jose Luis Mosquera}
%\note{}

\seealso{\code{\link{toOOC}}, }
\examples{
data(joslyn)

N.t <- Nt(joslyn.OOC)
print(N.t)
}
\keyword{ Nt }
\keyword{ toOOC }
