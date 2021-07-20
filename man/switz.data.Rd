\name{switz.data}
\alias{switz.data}
\docType{data}
\title{
Number of daily new cases of COVID-19 in Switzerland
}
\description{
A data set with the number of new cases of COVID-19 in Switzerland. A natural data set that conform very well to the Benford's law
}
\usage{data("switz.data")}
\format{
  A numeric vector with 383 observations.
  \describe{
    \item{\code{x}}{a numeric vector}
  }
}
\source{
Our World in Data COVID-19 project: https://ourworldindata.org/coronavirus-data
}
\references{
Roser M, Ritchie H, Ortiz-Ospina E, Hasel J. 2020. Coronavirus Pandemic (COVID-19). https://ourworldindata.org/coronavirus.}
\examples{
data(switz.data)
## maybe str(switz.data) ; plot(switz.data) ...
}
\keyword{datasets}
