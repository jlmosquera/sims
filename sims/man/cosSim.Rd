\name{cosSim}
\alias{cosSim}
\title{Cosine similarity measure}
\description{This function computes the cosine similarity measure between two columns of a \code{matrix}}
\usage{cosSim(x, na.rm = FALSE)}
\arguments{
  \item{x}{\code{matrix} with two columns}
  \item{na.rm}{a \code{logical} value indicating whether \code{NA} values should be stripped before the computation proceeds. By default is \code{FALSE}}
}
\details{Given the two columns of the \code{matrix}, let A and B be the names of such columns, the cosine similarity, is estimated using a dot product and magnitude as

\deqn{sim_{cos}(A,B) = \frac{A \cdot B}{\Vert A \Vert \cdot \Vert B \Vert} = \frac{\sum_{i=1}^{n}A_{i} \times B_{i}}{\sqrt{ \sum_{i=1}^{n} (A_{i})^{2}} \times \sqrt{ \sum_{i=1}^{n} (B_{i})^{2}}}}
}
%\value{}
%\references{}
\author{Jose Luis Mosquera}
%\note{}
%\seealso{\code{\link{}}}
\examples{
a <- sample(c(1:100, NA), 20)
idx.na <- sample(1:20, 4)
a[idx.na] <- NA
b <- sample(1:100, 20)
x <- as.matrix(cbind(a,b))

cos.sim <- cosSim(x, na.rm = TRUE)
print(cos.sim)

## Cosine similarity between semantic similarites of Resnik

data(prostateIds)                # Data set from the package goProfiles

pckg <- "org.Hs.eg.db"           # Organism package of humans

eg.we <- welsh01EntrezIDs[1:10]   # Entrez Gene IDs signature Welsh 01
eg.sg <- singh01EntrezIDs[1:10]   # Entrez Gene IDs signature Singh 01

WEvsSG.nb <- gosimsAvsB(eg1 = eg.we, eg2 = eg.sg, ontology = "MF", pckg = pckg, type = "nb", method = "Res")

idx.Inf <- which(WEvsSG.nb == Inf)
WEvsSG.nb[idx.Inf] <-NA

WEvsSG.cos.sim <- cosSim(WEvsSG.nb , na.rm = TRUE)
print(WEvsSG.cos.sim)
}
\keyword{ cosSim }
\keyword{ gosimsAvsB }
