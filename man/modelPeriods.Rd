\name{Spectral Periods for NGA Models}
\alias{modelPeriods}
\alias{getPeriod}
\title{Spectral Periods for NGA Models}
\description{
  The function \code{modelPeriods} returns a vector of periods (sec) for
  which the model coefficients are defined for the different NGA
  models.

  The function \code{getPeriod} determines whether or not a
  given period \code{T} has defined coefficients.  If not, the function
  returns the next-highest and next-lowest periods with defined
  coefficients.
}
\usage{
modelPeriods(model, positive = FALSE)
getPeriod(T, model)
}
\arguments{
  \item{model}{a string indicating the name of the model from which
    periods should be returned, i.e., \code{"AS08", "BA08", "CB08",} or
    \code{"CY08"}.}
  \item{positive}{logical value (\code{TRUE} or \code{FALSE}) indicating
    whether or not to return only positive periods (i.e., spectral
    periods 0.01 sec and greater), excluding PGA (\code{T = 0}), PGV
    (\code{T = -1}), and PGD (\code{T = -2}, in the case of CB08) from the
    list.}
  \item{T}{spectral period at which the ground motion calculation is to
    be performed (sec)} 
}
\details{
  The \code{modelPeriods} function is a generalization of
  \code{\link{periods.as},} \code{\link{periods.ba},}
  \code{\link{periods.cb},} and \code{\link{periods.cy}.}  The purpose
  of the \code{positive} argument is to separate spectral
  acceleration (which is a continuous function of \code{T}) from PGA,
  PGV, and PGD.  This is useful for interpolation purposes (only Sa may
  be interpolated) and for plotting the predicted response spectra.
}
\value{
  \code{modelPeriods} returns a vector of periods that have defined
  coefficients for the specified NGA model (sec).
  
  
  \code{getPeriods} returns a three-element list with components
  \code{interp,} \code{lower,} and \code{upper}:
  \tabular{lll}{
       \tab  \code{interp} \tab a logical value indicating whether or not
       interpolation is necessary given the spectral period \code{T} \cr
       \tab  \code{lower} \tab gives the greatest period less than
       \code{T} that has defined model coefficients (if \code{interp = TRUE}) \cr
       \tab  \code{upper} \tab gives the smallest period greater than
       \code{T} that has defined model coefficients (if \code{interp = TRUE})
  }
  If \code{interp = FALSE}, then \code{lower} and \code{upper} contain
  null values.
}
\seealso{
  \code{\link{periods}}, \code{\link{interpolate}}, \code{\link{Sa}},
  \code{\link{Sa.nga}}.
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
\examples{
# Example 1: List of periods for the AS08 model

# Entire list of ground motion parameters
modelPeriods(model = "AS08")

# List of spectral periods excluding PGA and PGV
modelPeriods(model = "AS08", positive = TRUE)


# Example 2: Find whether interpolation is necessary to estimate
# ground motions at a spectral period of 0.65 sec using the AS08 model

getPeriod(T = 0.65, model = "AS08")

# ANSWER: Yes, interpolation is necessary.  The next-lowest period
# with defined coefficients is 0.5 sec, and the next-highest period
# with defined coefficients is 0.75 sec.
}
