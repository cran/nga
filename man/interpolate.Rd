\name{Interpolation for Intermediate Spectral Periods}
\alias{interpolate}
\title{Interpolation for Intermediate Spectral Periods}
\description{
  Performs linear interpolation; simple wrapper of the internal R function
  \code{approx}.
}
\usage{
interpolate(x, x1, x2, y1, y2) 
}
\arguments{
  \item{x}{x-value at which the interpolated y-value is desired.} 
  \item{x1, x2}{two x-values.}
  \item{y1, y2}{two y-values.}
}
\details{
  For log-log interpolation, the arguments should be entered in log space.
}
\value{
  Interpolated value \code{y} corresponding to \code{x},
  using linear interpolation between points \code{(x1, y1)} and
  \code{(x2, y2)}.
}
\seealso{
  \code{\link{getPeriod}}, \code{\link{Sa}}, \code{\link{Sa.nga}}.
}
\author{James Kaklamanos <james.kaklamanos@tufts.edu> and
  Eric M. Thompson <eric.thompson@tufts.edu>}
\examples{
# Example interpolation of spectral acceleration:

# Assumed earthquake parameters:
M <- 6
Rjb <- 30
Vs30 <- 500
rake <- 90
epsilon <- 0

# Desired: Median Sa at T = 0.19 sec using the BA08 model

# Since there are no defined coefficients at T = 0.19 sec,
# log-log interpolation is necessary.

# First, find the periods directly above and below T = 0.19 sec
T1 <- getPeriod(T = 0.19, model = "BA08")$lower
T2 <- getPeriod(T = 0.19, model = "BA08")$upper
T1
T2

# Second, find the spectral accelerations for those periods
Sa1 <- Sa.ba(M = M, Rjb = Rjb, Vs30 = Vs30, rake = rake, epsilon =
epsilon, T = T1)
Sa2 <- Sa.ba(M = M, Rjb = Rjb, Vs30 = Vs30, rake = rake, epsilon =
epsilon, T = T2)
Sa1
Sa2

# Third, use the interpolate function to find Sa at T = 0.19 sec
# Note the use of log-log interpolation
LnSa <- interpolate(x = log(0.19), x1 = log(T1), x2 = log(T2),
                    y1 = log(Sa1), y2 = log(Sa2))
Sa <- exp(LnSa)
Sa
}
