\name{Estimation of Depth to Top of Rupture, Ztor}
\alias{Ztor.calc}
\title{Estimation of Depth to Top of Rupture, Ztor}
\description{
  Estimates the depth to top of rupture, Ztor.
}
\usage{
Ztor.calc(Zhyp = NA, W, dip, M = NA, rake = NA)
}
\arguments{
  \item{Zhyp}{Hypocentral depth of the earthquake (km)}
  \item{W}{Down-dip rupture width (km).}
  \item{dip}{Fault dip angle (deg).}
  \item{M}{Moment magnitude of earthquake.}
  \item{rake}{Rake angle of fault movement (deg).}
}
\details{
  To implement this function, W and dip must be specified, as well as
  either (1) Zhyp, or (2) M and rake, which are used to estimate Zhyp if
  Zhyp is unspecified (using correlations presented in Scherbaum et
  al. (2004).  Estimates of W and dip may be obtained using the
  functions \code{\link{W.calc}} and \code{\link{dip.calc}},
  respectively.  The estimate of the depth to top of rupture (Ztor) is
  then obtained from Zhyp, W, and dip.  This calculation assumes that
  the hypocenter is located 60\% down the fault width, as suggested by
  Mai et al. (2005).
  
  Another currently unimplemented option is to estimate Ztor directly
  using a simple correlation between Ztor and M.  The depth to top of
  rupture is estimated using a quadratic relationship derived from the
  correlation suggested by Campbell et al. (2009).
}
\value{
  Estimated value of the depth to top of rupture, Ztor (km).
}
\seealso{
  \code{\link{W.calc}}, \code{\link{dip.calc}}, \code{\link{Sa}},
  \code{\link{Sa.nga}}.
}
\author{James Kaklamanos <james.kaklamanos@tufts.edu> and
  Eric M. Thompson <eric.thompson@tufts.edu>}
\references{
  
  Campbell, K., N. Abrahamson, M. Power, B. Chiou, Y. Bozorgnia,
  T. Shantz, and C. Roblee (2009).  Next Generation Attenuation (NGA)
  Project: Empirical Ground Motion Prediction Equations for Active
  Tectonic Regions. \emph{Sixth International Conference on Urban
    Earthquake Engineering}, Tokyo, Japan.
  
  Mai, P. M., P. Spudich, and J. Boatwright (2005).  Hypocenter
  Locations in Finite-Source Rupture Models. \emph{Bulletin of the
    Seismological Society of America} \bold{95,} 965--980.
}
\examples{
# Assumed earthquake parameters for this example:
Zhyp <- 8
M <- 6
rake <- 180   # Strike-slip fault

# Example 1:  Estimation of Ztor from Zhyp, W, and dip

# First, estimate W using W.calc
W <- W.calc(M = M, rake = rake)
W

# Second, estimate dip using dip.calc
dip <- dip.calc(rake = rake)
dip

# Third, estimate Ztor (now that we have estimates of W and dip)
Ztor1 <- Ztor.calc(Zhyp = Zhyp, W = W, dip = dip)
Ztor1


# Example 2:  Estimation of Ztor when Zhyp is not specified

# Say that we do not have an estimate for Zhyp, and that we wish to
# estimate Ztor.  In lieu of Zhyp, we use M and rake as inputs
# to the function.

Ztor2 <- Ztor.calc(W = W, dip = dip, M = M, rake = rake)
Ztor2
}
