\name{Estimation of Fault Dip}
\alias{dip.calc}
\title{Estimation of Fault Dip}
\description{
  Estimates the fault dip angle from the style of faulting (using the
  rake angle), following the explanation in Kaklamanos et al. (2011).
  These recommendations are a modification of the guidelines Chiou and
  Youngs (2008) utilized in developing their NGA model.
}
\usage{
dip.calc(rake)
}
\arguments{
  \item{rake}{Rake angle of fault movement (deg).}
}
\value{
  Estimated fault dip angle (deg).
}
\seealso{
  \code{\link{Ztor.calc}}, \code{\link{W.calc}}, \code{\link{Sa}},
  \code{\link{Sa.nga}}.
}
\author{James Kaklamanos <james.kaklamanos@tufts.edu> and
  Eric M. Thompson <eric.thompson@tufts.edu>}
\references{
  Chiou, B. S.-J., and R. R. Youngs (2008).  NGA Model for the Average
  Horizontal Component of Peak Ground Motion and Response Spectra,
  \emph{PEER Report No. 2008/09}, Pacific Earthquake Engineering
  Research Center, University of California, Berkeley.
  
  Kaklamanos, J., L. G. Baise, and D. M. Boore (2011).  Estimating
  Unknown Input Parameters when Implementing the NGA Ground-Motion
  Prediction Equations in Engineering Practice.
  \emph{Earthquake Spectra} \bold{27,} 1219--1235.
}
\examples{
# Estimated dip angle for a strike-slip fault
dip.calc(rake = 180)

# Estimated dip angle for a reverse fault
dip.calc(rake = 90)

# Estimated dip angle for a normal fault
dip.calc(rake = -90) 
}
