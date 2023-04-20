\name{ancestors}
\alias{ancestors}
\title{Ancestors for each term of the ontology}
\description{Given the accesibility matrix, generates a list whose elements are the ancestors for each term of the ontology}
\usage{ancestors(m)}
\arguments{
	\item{m}{\code{matrix} of accessibility}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{The resulting object is a \code{list} whose elements are \code{character} 
vectors with the names of the ancestors associated with each term
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{Jose Luis Mosquera}
%\note{
%%  ~~further notes~~
%}

\seealso{\code{\link{getA}}}
\examples{
data(joslyn)

inv.IminusG <- inverseIminusG(joslyn.OOC)
A.mat <- getA(inv.IminusG)

at <- ancestors(A.mat)
print(at)
}
\keyword{ ancestors }
\keyword{ getA }
\keyword{ accesibility matrix}
