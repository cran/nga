\name{AS08 Sub-Functions}
\alias{subs.as}
\alias{f1.as}
\alias{f5.as}
\alias{f4.as}
\alias{f6.as}
\alias{f8.as}
\alias{f10.as}
\alias{Alpha.as}
\alias{Sigma0.as}
\alias{Tau0.as}
\alias{Sigma.as}
\alias{Tau.as}
\alias{SigmaTot.as}
\alias{PGA1100.as}
\alias{SaMedian.as}
\title{AS08 Sub-Functions}
\description{
  The sub-functions needed to estimate ground motions using the
  Abrahamson & Silva NGA model.  These functions are called by the main
  function used for calculations, \code{\link{Sa.as}}.
}
\usage{
f1.as(M, Rrup, T)
f5.as(PGA1100, Vs30, T)
f4.as(Rjb, Rx, dip, Ztor, M, W, T)
f6.as(Ztor, T)
f8.as(Rrup, M, T)
f10.as(Z1.0, Vs30, T)
PGA1100.as(M, Rrup, Rjb, Rx, Ztor, Frv, Fnm, Fas, Fhw, dip, W)
Alpha.as(PGA1100, Vs30, T)
Sigma0.as(M, VsFlag, T)
Tau0.as(M, T)
Sigma.as(M, PGA1100, Vs30, VsFlag, T)
Tau.as(M, PGA1100, Vs30, T)
SigmaTot.as(M, PGA1100, Vs30, VsFlag, T)
SaMedian.as(M, Rjb, Rrup, Rx, Ztor, Frv, Fnm, Fas, Fhw, dip, Vs30, Z1.0, W, T)
}
\arguments{
  \item{M}{Moment magnitude of earthquake.}
  \item{Rjb}{Horizontal distance to the surface projection of the
    rupture plane; Joyner-Boore distance (km).}
  \item{Rrup}{Closest distance to the rupture plane; rupture distance
    (km).}
  \item{Rx}{Horizontal distance to the surface projection of the
    top edge of the rupture plane, measured perpendicular to the
    strike; site coordinate (km).}
  \item{Ztor}{Depth to top of rupture (km).}
  \item{Frv}{Reverse style-of-faulting flag; equal to 1 for reverse
    faulting (\code{30 <= rake <= 150 deg}), and 0 otherwise.}
  \item{Fnm}{Normal style-of-faulting flag; equal to 1 for normal
    faulting (\code{-120 <= rake <= -60 deg}), and 0 otherwise.}
  \item{Fas}{Aftershock flag; equal to 1 for aftershocks and 0 for
    mainshocks.}
  \item{Fhw}{Hanging wall flag; equal to 1 for sites on the hanging
    wall side of the fault (\code{Rx >= 0} and \code{dip != 90}), and 0
    otherwise.}
  \item{dip}{Fault dip angle (deg).}
  \item{Vs30}{Time-averaged shear wave velocity over a subsurface
    depth of 30 meters (m/s).}
  \item{VsFlag}{Flag variable indicating how \code{Vs30} is obtained;
    equal to 1 if \code{Vs30} is measured, and 0 if \code{Vs30} is
    estimated or inferred.}
  \item{Z1.0}{Depth to Vs = 1.0 km/s (m).}
  \item{W}{Down-dip rupture width (km).}
  \item{PGA1100}{Median PGA when \code{Vs30 = 1100 m/s} (g).}
  \item{T}{Spectral period (sec). Use 0 for PGA and -1 for PGV.}
}
\details{
  \describe{
    \item{\code{f1.as}}{Base model (Eqns 2 and 3 of Abrahamson and
    Silva (2008)).}
    \item{\code{f5.as}}{Site response model (Eqns 4 to 6).}
    \item{\code{f4.as}}{Hanging wall model (Eqns 7 to 12); Eqn 12
    modified in Abrahamson and Silva (2009).}
    \item{\code{f6.as}}{Depth to top of rupture model (Eqn 13).}
    \item{\code{f8.as}}{Large distance model (Eqns 14 and 15).}
    \item{\code{f10.as}}{Soil depth model (Eqns 16 to 20).}
    \item{\code{Alpha.as}}{Partial derivative of site response
    function (\code{f5.as}) with respect to \code{PGA1100} (Eqn 26;
    corrected in Abrahamson and Silva (2009)).}
    \item{\code{Sigma0.as}}{Intra-event standard deviation of the
    observed ground motion for the linear site response range (Eqn
    27).}
    \item{\code{Tau0.as}}{Inter-event standard deviation of the observed
    ground motion for the linear site response range (Eqn 28).}
    \item{\code{Sigma.as}}{Intra-event standard deviation (Eqn 24;
    corrected in Abrahamson and Silva (2009)).}
    \item{\code{Tau.as}}{Inter-event standard deviation (Eqn 25).}
    \item{\code{SigmaTot.as}}{Total standard deviation.}
    \item{\code{PGA1100.as}}{Calculation of median PGA when \code{Vs30
    = 1100 m/s}}
    \item{\code{SaMedian.as}}{Median estimate of ground motion (Eqn 1).}
  }
}
\value{
  Returns values of the sub-functions for the Abrahamson & Silva
  NGA model.  The outputs of the functions are in log space, with
  the exception of (1) \code{SaMedian.as}, which outputs values in g
  (for Sa and PGA) or cm/s (for PGV), and (2) \code{PGA1100.as}, which
  outputs values in g.
}
\seealso{
  \code{\link{coefs}}, \code{\link{coefs.t.as}},
  \code{\link{periods.as}}, \code{\link{Sa.as}}, \code{\link{Sa.nga}}.
}
\author{James Kaklamanos <james.kaklamanos@tufts.edu> and
  Eric M. Thompson <eric.thompson@tufts.edu>}
\references{
  Abrahamson, N., and W. Silva (2008). Summary of the Abrahamson & Silva
  NGA Ground-Motion Relations. \emph{Earthquake Spectra} \bold{24,}
  67--97.
   
  Abrahamson, N., and W. Silva (2009). Errata for \dQuote{Summary
    of the Abrahamson & Silva NGA Ground-Motion Relations} by
  Abrahamson, N.A. and W.J. Silva (2008). Pacific Earthquake Engineering
  Research Center, University of California, Berkeley,
  <http://peer.berkeley.edu/products/nga_project.html>.
}
\keyword{internal}
