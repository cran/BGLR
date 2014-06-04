\name{effects.BGLR}
\alias{effects.BGLR}
\title{effects for BGLR fitted models}
\description{
    Gives estimated effects for models fitted using the BGLR function.
}
\usage{
\method{effects}{BGLR}(object, ...)
}
\arguments{
  \item{object}{An object of class \code{BGLR}.}
  \item{...}{Further arguments passed to or from other methods.}
}
\author{
  Gustavo de los Campos, Paulino Perez Rodriguez,
}
\seealso{
  \code{BGLR}.
}
\examples{

\dontrun{

setwd(tempdir())
library(BGLR)
data(wheat)
out=BLR(y=wheat.Y[,1],XL=wheat.X)
effects(out)

}

}
\keyword{summary}
