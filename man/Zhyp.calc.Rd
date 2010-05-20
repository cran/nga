\name{Estimation of Hypocentral Depth, Zhyp}
\alias{Zhyp.calc}
\title{Estimation of Hypocentral Depth, Zhyp}
\description{
  Provides an estimate of the hypocentral depth (\code{Zhyp}), which in
  turn may be used to estimate the depth to top of rupture
  (\code{Ztor}), if \code{Ztor} is unknown.
}
\usage{
Zhyp.calc(M, rake)
}
\arguments{
  \item{M}{Moment magnitude of earthquake.}
  \item{rake}{Rake angle of fault movement (deg).}
}
\details{
  The value of \code{Zhyp} is estimated using correlations presented in
  Table 1 of Scherbaum et al. (2004).
}
\value{
  Estimated value of the hypocentral depth, \code{Zhyp} (km).
}
\seealso{
  \code{\link{Ztor.calc}}, \code{\link{Sa}}, \code{\link{Sa.nga}}.
}
\author{James Kaklamanos <james.kaklamanos@tufts.edu> and
  Eric M. Thompson <eric.thompson@tufts.edu>}
\references{
  Scherbaum, F., J. Schmedes, and F. Cotton (2004).  On the Conversion
  of Source-to-Site Distance Measures for Extended Earthquake Source
  Models. \emph{Bulletin of the Seismological Society of America}
  \bold{94,} 1053--1069.
}
\examples{
# Estimate the hypocentral depths for some various scenarios:

# Small earthquake, shallow-dipping fault (not strike-slip)
Zhyp.calc(M = 5, rake = 90)

# Small earthquake, strike-slip fault
Zhyp.calc(M = 5, rake = 180)

# Large earthquake, shallow-dipping fault (not strike-slip)
Zhyp.calc(M = 7, rake = -90)

# Large earthquake, strike-slip fault
Zhyp.calc(M = 7, rake = 0)
}