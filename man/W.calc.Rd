\name{Estimation of Down-Dip Rupture Width, W}
\alias{W.calc}
\title{Estimation of Down-Dip Rupture Width, W}
\description{
  Estimates the down-dip rupture width (\code{W}) from the moment
  magnitude of the earthquake (\code{M}) using the empirical
  correlations published in Wells and Coppersmith (1994), for
  strike-slip, normal, and reverse faulting mechanisms.
}
\usage{
W.calc(M, rake)
}
\arguments{
  \item{M}{Moment magnitude of the earthquake.}
  \item{rake}{Rake angle of fault movement (deg).}
}
\value{
  Estimated down-dip width of the rupture plane, \code{W} (km).
}
\seealso{
  \code{\link{Ztor.calc}}, \code{\link{dip.calc}}, \code{\link{Sa}},
  \code{\link{Sa.nga}}.
}
\author{James Kaklamanos <james.kaklamanos@tufts.edu> and
  Eric M. Thompson <eric.thompson@tufts.edu>}
\references{
  Wells, D. L., and K. J. Coppersmith (1994).  New Empirical
  Relationships among Magnitude, Rupture Length, Rupture Width, Rupture
  Area, and Surface Displacement. \emph{Bulletin of the Seismological
    Society of America} \bold{84,} 974--1002.
}
\examples{
# Estimate the down-dip rupture widths for some various scenarios

# Small earthquake, reverse fault
W.calc(M = 5, rake = 90)

# Small earthquake, normal fault
W.calc(M = 5, rake = -90)

# Small earthquake, strike-slip fault
W.calc(M = 5, rake = 180)

# Large earthquake, reverse fault
W.calc(M = 7, rake = 90)

# Large earthquake, strike-slip fault
W.calc(M = 7, rake = 0)

# Large earthquake, normal fault
W.calc(M = 7, rake = 90)
}
