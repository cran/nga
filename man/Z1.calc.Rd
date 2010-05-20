\name{Estimation of Depth Parameter, Z1.0}
\alias{Z1.calc}
\alias{Z1.calc.as}
\alias{Z1.calc.cy}
\title{Estimation of Depth Parameter, Z1.0}
\description{
  Estimates the depth parameter \code{Z1.0} from the average shear wave
  velocity (\code{Vs30}), using Equation 17 in Abrahamson and Silva
  (2008) and Equation 1 in Chiou and Youngs (2008) for \code{Z1.calc.as}
  and \code{Z1.calc.cy}, respectively.
}
\usage{
Z1.calc.as(Vs30)
Z1.calc.cy(Vs30)
}
\arguments{
  \item{Vs30}{Time-averaged shear wave velocity over a subsurface
    depth of 30 meters (m/s).}
}
\value{
  Estimated value of \code{Z1.0}, the depth to a shear wave velocity
  horizon of Vs = 1.0 km/s  (m).
}
\seealso{
  \code{\link{Sa}}, \code{\link{Sa.nga}}.
}
\author{James Kaklamanos <james.kaklamanos@tufts.edu> and
  Eric M. Thompson <eric.thompson@tufts.edu>}
\references{
  Abrahamson, N., and W. Silva (2008). Summary of the Abrahamson & Silva
  NGA Ground-Motion Relations. \emph{Earthquake Spectra} \bold{24,}
  67--97.

  Chiou, B. S.-J., and R. R. Youngs (2008). An NGA Model for the Average
  Horizontal Component of Peak Ground Motion and Response Spectra.
  \emph{Earthquake Spectra} \bold{24,} 173--215.
}
\examples{
# Estimated depth to Vs = 1.0 km/s using the AS08 and CY08 correlations

# AS08 model, Vs30 = 500 m/s
Z1.calc.as(Vs30 = 500)

# CY08 model, Vs30 = 500 m/s
Z1.calc.cy(Vs30 = 500)

# The CY08 relation generates smaller values of Z1.0 than the
# AS08 relation generates.
}
