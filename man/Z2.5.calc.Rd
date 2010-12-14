\name{Estimation of Depth Parameter, Z2.5}
\alias{Z2.5.calc}
\title{Estimation of Depth Parameter, Z2.5}
\description{
  Estimates the depth parameter \code{Z2.5} from either \code{Z1.5},
  \code{Z1.0}, or \code{Vs30}, in decreasing order of preference, using
  the guidelines by Campbell and Bozorgina (2007)
}
\usage{
Z2.5.calc(Vs30 = NA, Z1.0 = NA, Z1.5 = NA)
}
\arguments{
  \item{Vs30}{Time-averaged shear wave velocity over a subsurface
    depth of 30 meters (m/s).}
  \item{Z1.0}{Depth to Vs = 1.0 km/s (m).}
  \item{Z1.5}{Depth to Vs = 1.5 km/s (m).}
}
\value{
  Estimated value of \code{Z2.5}, the depth to a shear wave velocity
  horizon of Vs = 2.5 km/s  (m).
}
\seealso{
  \code{\link{Sa}}, \code{\link{Sa.nga}}, \code{\link{Z1.calc}}.
}
\author{James Kaklamanos <james.kaklamanos@tufts.edu> and
  Eric M. Thompson <eric.thompson@tufts.edu>}
\references{
  Campbell, K. W., and Y. Bozorgnia (2007). Campbell-Bozorgnia NGA
  Ground Motion Relations for the Geometric Mean Horizontal Component of
  Peak and Spectral Ground Motion Parameters, \emph{PEER Report
    No. 2007/02}, Pacific Earthquake Engineering Research Center,
  University of California, Berkeley.
}
\examples{
# Estimated depth to Vs = 2.5 km/s

# Example if Z1.5 is known
Z2.5.calc(Z1.5 = 1000)

# Example if Z1.0 is known
Z2.5.calc(Z1.0 = 800)

# Example if only Vs30 is known
Z2.5.calc(Vs30 = 400)

}
