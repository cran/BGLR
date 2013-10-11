\name{predict.BGLR}
\alias{predict.BGLR}
\title{Predictions from BGLR Analysis}
\description{
    Predicting values using results from BGLR function.
}
\usage{
\method{predict}{BGLR}(object,newdata = NULL, ...)
}
\arguments{
  \item{object}{An object of class \code{BGLR}.}
  \item{newdata}{new data, see \code{BGLR} function for more details.}
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

}

}
\keyword{regression}
