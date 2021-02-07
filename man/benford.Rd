\name{benford}
\alias{benford}
\title{Benford's analysis
}
\description{
Benford's analysis makes use of a statistic property of natural data sets called Benford's Law. Benford’s Law (also called “first digit phenomenon”) is a statistical phenomenon that describes the frequency of a given integer, from 1 to 9, to be in the first significant digit in the numbers of a large data set. The Benford’s law has been most practically used to detect fraud or rounding errors in real world numbers. This is possible by examining departures in the frequencies of individual digits from those predicted by Benford. This only makes sense once it is established (often empirically) that
the data follow the law under normal circumstances (Sambridge et al., 2011). This is true because human pseudo-random productions are in many ways different from true randomness (Nickerson, 2002). As a consequence, fabricated data might fit to the Benford’s Law to a lesser extent than genuine data (Banks and Hill, 1974; Gauvrit et al., 2017).
}
\usage{
benford(x, plot = FALSE)
}
\arguments{
  \item{x}{
A numeric vector with the data set numbers to be analyzed
}
  \item{plot}{
A logic that control whether the resulting first digit distribution and the Benford's distribution would be ploted
}
}
\value{
LIST countaining:
1. Named vector with three elements: the Chi Square test p value (p), the root mean square deviation (RMSD) from the Benford's distribution, and the log of the likelihhod of the first digit distribution in relation to the Benford's distribution;
2. Matrix, with three columns, countaining the first digits ([,1]), the frequency counts of the first digit in the data set ([,2]) and the frequency count of the first digit in a classic Benford's distribution ([,3])
}
\references{
Banks WP, Hill DK. 1974. The apparent magnitude of number scaled by random production. J. Exp. Psychol. 102:353–376. http://content.apa.org/journals/xge/102/2/353.
Benford F. 1938. The Law of Anomalous Numbers. Proc. Am. Philos. Soc. 78:551–572. http://www.jstor.org/stable/984802.
Gauvrit N, Houillon J-C, Delahaye J-P. 2017. Generalized Benford’s Law as a Lie Detector. Adv. Cogn. Psychol. 13:121–127. http://ac-psych.org/en/download-pdf/id/214.
Nickerson RS. 2002. The production and perception of randomness. Psychol. Rev. 109:330–357. http://doi.apa.org/getdoi.cfm?doi=10.1037/0033-295X.109.2.330.
Sambridge M, Tkalcic H, Arroucau P. 2011. Benford’s Law of First Digits: From Mathematical Curiosity to Change Detector. Asia Pacific Math. Newsl. 1:1–5.
}
\author{
Vitor Hugo Moreau, Ph.D.
Department of Biotechnology
Federal University of Bahia, Brazil
}
\note{
RMSD and likelihood are not formal statistic tests, so it may be evaluated only in a comparative way. To perform analysis in order to get to absolute conclusion on the veracity of the data set, Chi square p value is more trustable.
}
\examples{
#Computer generated random data do not conform to the benford law
result <- benford(seq(1,10000)+rnorm(10000,0,100), FALSE)
#Natural data set, countaining the number of daily new cases of COVID-19 in Switzerland
##conform to the Benford' Law
result <- benford(switz.data, TRUE)
}
\keyword{ Benford's analysis }
\keyword{ Benford's Law }
