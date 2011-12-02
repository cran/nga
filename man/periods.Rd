\name{Spectral Periods for Each NGA Model}
\alias{periods}
\alias{periods.as}
\alias{periods.ba}
\alias{periods.cb}
\alias{periods.cy}
\title{Spectral Periods for Each NGA Model}
\description{
  Returns a vector of periods (sec) at which the period-dependent model
  coefficients are defined for the different NGA models.
}
\usage{
periods.as(positive = FALSE)
periods.ba(positive = FALSE)
periods.cb(positive = FALSE)
periods.cy(positive = FALSE)
}
\arguments{
  \item{positive}{logical value (\code{TRUE} or \code{FALSE}) indicating
    whether or not to return only positive periods (i.e., spectral
    periods 0.01 sec and greater), excluding PGA (\code{T = 0}), PGV
    (\code{T = -1}), and PGD (\code{T = -2}, in the case of CB08) from the
    list.  If \code{positive = FALSE}, the periods corresponding to
    PGA, PGV, and PGD (for the CB08 model) are appended to the list.}
}
\details{
  The purpose of the \code{positive} argument is to separate spectral
  acceleration (which is a continuous function of \code{T}) from PGA,
  PGV, and PGD.  This is useful for interpolation purposes (only Sa may
  be interpolated) and for plotting the predicted response spectra.
  \tabular{ll}{
    \code{periods.as} \tab Returns an array of periods for the AS08
    model. \cr
    \code{periods.ba} \tab Returns an array of periods for the BA08
    model. \cr
    \code{periods.cb} \tab Returns an array of periods for the CB08
    model. \cr
    \code{periods.cy} \tab Returns an array of periods for the CY08
    model. \cr
  }
}
\value{
  Vector of periods that have defined model coefficients (sec). 
}
\seealso{
  \code{\link{coefs}}, \code{\link{coefs.t.as}},
  \code{\link{coefs.t.ba}}, \code{\link{coefs.t.cb}},
  \code{\link{coefs.t.cy}}, \code{\link{modelPeriods}},
  \code{\link{getPeriod}}.
}
\author{James Kaklamanos <james.kaklamanos@tufts.edu> and
  Eric M. Thompson <eric.thompson@tufts.edu>}
\references{
  Abrahamson, N., and W. Silva (2008). Summary of the Abrahamson & Silva
  NGA Ground-Motion Relations. \emph{Earthquake Spectra} \bold{24,}
  67--97.
  
  Boore, D. M., and G. M. Atkinson (2008). Ground-Motion Prediction
  Equations for the Average Horizontal Component of PGA, PGV, and
  5\%-Damped PSA at Spectral Periods between 0.01 s and 10.0
  s. \emph{Earthquake Spectra} \bold{24,} 99--138.
  
  Campbell, K. W., and Y. Bozorgnia (2008). NGA Ground Motion Model for
  the Geometric Mean Horizontal Component of PGA, PGV, PGD, and 5\%
  Damped Linear Elastic Response Spectra for Periods Ranging from 0.01
  to 10 s. \emph{Earthquake Spectra} \bold{24,} 139--171.

  Chiou, B. S.-J., and R. R. Youngs (2008). An NGA Model for the Average
  Horizontal Component of Peak Ground Motion and Response Spectra.
  \emph{Earthquake Spectra} \bold{24,} 173--215.
}
\keyword{internal}
