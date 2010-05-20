\name{Estimation of Depth to Top of Rupture, Ztor}
\alias{Ztor.calc}
\title{Estimation of Depth to Top of Rupture, Ztor}
\description{
  Estimates the depth to top of rupture, \code{Ztor}.
}
\usage{
Ztor.calc(W, dip, Zhyp)
}
\arguments{
  \item{W}{Down-dip rupture width (km).}
  \item{dip}{Fault dip angle (deg).}
  \item{Zhyp}{Hypocentral depth of the earthquake (km).}
}
\details{
  To implement this function, \code{W}, \code{dip}, and \code{Zhyp} must
  be specified.  Estimates of \code{W}, \code{dip}, and \code{Zhyp} may
  be obtained using the functions \code{\link{W.calc}},
  \code{\link{dip.calc}}, and \code{\link{Zhyp.calc}}, respectively.
  The resulting calculation for \code{Ztor} assumes that the hypocenter
  is located 60\% down the fault width, as suggested by Mai et al. (2005).
}
\value{
  Estimated value of the depth to top of rupture, Ztor (km).
}
\seealso{
  \code{\link{W.calc}}, \code{\link{dip.calc}}, \code{\link{Zhyp.calc}},
  \code{\link{Sa}}, \code{\link{Sa.nga}}.
}
\author{James Kaklamanos <james.kaklamanos@tufts.edu> and
  Eric M. Thompson <eric.thompson@tufts.edu>}
\references{
  Mai, P. M., P. Spudich, and J. Boatwright (2005).  Hypocenter
  Locations in Finite-Source Rupture Models. \emph{Bulletin of the
    Seismological Society of America} \bold{95,} 965--980.
}
\examples{
# Assumed earthquake parameters for this example:
M <- 6
rake <- 180   # Strike-slip fault

# First, estimate W using W.calc
W <- W.calc(M = M, rake = rake)
W

# Second, estimate dip using dip.calc
dip <- dip.calc(rake = rake)
dip

# Third, estimate Zhyp using Zhyp.calc
Zhyp <- Zhyp.calc(M = M, rake = rake)
Zhyp

# Third, estimate Ztor (now that we have estimates of W, dip, and Zhyp)
Ztor <- Ztor.calc(W = W, dip = dip, Zhyp = Zhyp)
Ztor
}