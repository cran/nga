\name{Ground Motion Predictions for Individual Models}
\alias{Sa.as}
\alias{Sa.ba}
\alias{Sa.cb}
\alias{Sa.cy}
\alias{Sa}
\title{Ground Motion Predictions for Individual Models}
\description{
  Main functions for estimating ground motion parameters using the
  ground motion prediction equations developed during the Next
  Generation Attenuation of Ground Motions (NGA) project in 2008.  Each
  function performs ground motion calculations using for an individual
  NGA model.
}
\usage{
Sa.as(M, Rjb, Vs30, VsFlag, epsilon, T, Rrup = NA, Rx = NA,
      dip = NA, W = NA, Ztor = NA, Z1.0 = NA, rake = NA, Frv = NA,
      Fnm = NA, Fhw = NA, azimuth = NA, Zhyp = NA, Fas = 0)
Sa.ba(M, Rjb, Vs30, epsilon, T, rake = NA, U = NA, SS = NA,
      NS = NA, RS = NA, AB11 = 0)
Sa.cb(M, Rjb, Vs30, epsilon, T, Rrup = NA, dip = NA, W = NA,
      Ztor = NA, Z1.0 = NA, Z1.5 = NA, Z2.5 = NA, rake = NA,
      Frv = NA, Fnm = NA, Fhw = NA, azimuth = NA, Zhyp = NA, arb = 0)
Sa.cy(M, Rjb, Vs30, VsFlag, epsilon, T, Rrup = NA, Rx = NA,
      dip = NA, W = NA, Ztor = NA, Z1.0 = NA, rake = NA, Frv = NA,
      Fnm = NA, Fhw = NA, azimuth = NA, Zhyp = NA, AS = 0)
}
\arguments{
  \item{M}{Moment magnitude of earthquake.}
  \item{Rjb}{Joyner-Boore distance (km):  the horizontal distance to the
    surface projection of the rupture plane.}
  \item{Vs30}{Time-averaged shear wave velocity over a subsurface
    depth of 30 meters (m/s).}
  \item{VsFlag}{Flag variable indicating how \code{Vs30} is obtained;
    equal to 1 if \code{Vs30} is measured, and 0 if \code{Vs30} is
    estimated or inferred.}
  \item{epsilon}{number of standard deviations to be considered in the
    calculations.  Use 0 to obtain a median estimate of ground motion.}
  \item{T}{Spectral period (sec).  Use 0 for PGA and -1 for PGV. For the
    CB08 model only, specify -2 for PGD.  For spectral acceleration,
    \code{T} must be in the range \code{0.01 <= T <= 10 sec}.  If the
    specified period is within the allowable range and does not have
    defined equations, the program uses log-log interpolation (using
    \code{\link{interpolate}}) between the next-highest and next-lowest
    spectral periods with defined equations.}
  \item{Rrup}{Rupture distance (km):  the closest distance to the
    rupture plane; if left empty, \code{Rrup} is calculated from
    \code{Rx}, the source-to-site azimuth, and the geometric rupture
    parameters (\code{Ztor}, \code{W}, and \code{dip}) using
    \code{\link{Rrup.calc}}.}
  \item{Rx}{Site coordinate (km):  The horizontal distance to the
    surface projection of the top edge of the rupture plane, measured
    perpendicular to the strike.  If left empty, \code{Rx} is calculated
    from \code{Rjb}, the source-to-site azimuth, and the geometric
    rupture parameters (\code{Ztor}, \code{W}, and \code{dip}) using
    \code{\link{Rx.calc}}.  When only \code{Rjb} and the azimuth are
    assumed, \code{Rjb} is used to calculate \code{Rx}, which is then
    used to calculate \code{Rrup}.}
  \item{dip}{Dip angle of the rupture plane (deg).  If left empty, the
    dip is estimated using \code{\link{dip.calc}}.}
  \item{W}{Down-dip width of rupture plane (km).  If left empty,
    \code{W} is estimated using \code{\link{W.calc}}.}
  \item{Ztor}{Depth to top of rupture (km).  If left empty, \code{Ztor}
    is estimated using \code{\link{Ztor.calc}}.}
  \item{Z1.0}{Depth to Vs = 1.0 km/s (m).  If left empty, \code{Z1.0} is
    estimated using \code{\link{Z1.calc.as}} for the AS08 model and
    \code{\link{Z1.calc.cy}} for the CY08 model.}
  \item{Z1.5}{Depth to Vs = 1.5 km/s (m).  Z1.5 is not utilized in
    ground motion calculations, but if available, it may be used to
    estimate \code{Z2.5} for the CB08 model.}
  \item{Z2.5}{Depth to Vs = 2.5 km/s (m; \bold{note the units}). If left
    empty, \code{Z2.5} is estimated from \code{Z1.5} or \code{Z1.0} if
    available, using the recommendations in Campbell and Bozorgnia
    (2007).  If neither \code{Z1.5} nor \code{Z1.0} is available, then
    \code{Vs30} is used to estimate \code{Z1.0} using
    \code{\link{Z1.calc.as}}, which is in turn used to estimate \code{Z2.5}.}
  \item{rake}{Rake angle of fault movement (deg).  Either the rake angle
    or the style-of-faulting flag variables (\code{Frv} and \code{Fnm}
    for AS08, CB08, and CY08; and \code{U}, \code{RS}, \code{NS}, and
    \code{SS} for BA08) must be specified.}
  \item{Frv}{Reverse style-of-faulting flag (1 for reverse faulting, 0
    otherwise) for the AS08, CB08, and CY08 models.  Either (a) the rake
    angle, or (b) both \code{Frv} and \code{Fnm}, must be specified.
    Reverse faulting is characterized by rake angles in the range
    \code{30 <= rake <= 150 deg} for the AS08 and CY08 models, and in
    the range \code{30 < rake < 150 deg} for the CB08 model.}
  \item{Fnm}{Normal style-of-faulting flag (1 for normal faulting, 0
    otherwise) for the AS08, CB08, and CY08 models.  Either (a) the rake
    angle, or (b) both \code{Frv} and \code{Fnm}, must be specified.
    Normal faulting is characterized by rake angles in the range
    \code{-120 <= rake <= -60 deg} for the AS08 and CY08 models and
    \code{-150 < rake < -30 deg} for the CB08 model.}
  \item{U}{Unspecified style-of-faulting flag for the BA08 model;
    equal to 1 if the user wishes to perform a generic ground motion
    calculation when the style of faulting is unspecified, and 0
    otherwise.}
  \item{RS}{Reverse style-of-faulting flag for the BA08 model;
    equal to 1 for reverse faulting (\code{30 <= rake <= 150 deg}), and
    0 otherwise.}
  \item{NS}{Normal style-of-faulting flag for the BA08 model;
    equal to 1 for normal faulting (\code{-150 <= rake <= -30 deg}), and
    0 otherwise.}
  \item{SS}{Strike-slip style-of-faulting flag for the BA08 model; equal
    to 1 for strike-slip faulting (when the rake is not in either of the
    ranges specified for RS or NS), and 0 otherwise.}
  \item{Fhw}{Hanging wall flag; equal to 1 for sites on the hanging
    wall side of the fault (\code{Rx >= 0; azimuth >= 0}), and 0
    otherwise.  For AS08 and CY08, either \code{Fhw}, \code{Rx}, or
    \code{azimuth} must be specified.  For CB08, the parameters
    \code{Fhw} and \code{azimuth} are optional, and they are only used
    to estimate \code{Rrup} when \code{Rrup} is unknown; if neither
    \code{Fhw} nor \code{azimuth} is specified, the site is assumed to
    be located on the footwall, and \code{Rrup} is easily estimated as
    \code{sqrt(Rjb^2 + Ztor^2).}}
  \item{azimuth}{Source-to-site azimuth (deg); see Kaklamanos et
    al. (2011) for details.  Used by \code{\link{Rx.calc}} and
    \code{\link{Rrup.calc}} for distance calculations.}
  \item{Zhyp}{Hypocentral depth of the earthquake (km).  \code{Zhyp} is not
    utilized in ground motion calculations, but it may be used to estimate
    \code{Ztor}.  See \code{\link{Ztor.calc}} for details.}
  \item{Fas}{Aftershock flag for AS08; equal to 1 for aftershocks and 0 for
    mainshocks (the default).}
  \item{AS}{Aftershock flag for CY08; equal to 1 for aftershocks and 0 for
    mainshocks (the default).}
  \item{arb}{Flag variable indicating the method of determining aleatory
    uncertainty for the CB08 model; equal to 1 if the standard deviation
    should be calculated for the arbitrary horizontal component of ground
    motion, and 0 if the standard deviation should be calculated for the
    geometric mean horizontal ground motion (the default).}
  \item{AB11}{Flag variable equaling 1 if the Atkinson and Boore (2011)
    small-magnitude correction factor should be applied to the BA08
    model, and 0 otherwise.}
  
}
\value{
  The spectral acceleration (in units of \emph{g}) at period T; peak
  ground acceleration (PGA, in units of \emph{g}) when \code{T = 0};
  peak ground velocity (PGV, in units of cm/sec) when \code{T = -1}; and
  peak ground displacement using the CB08 model (PGD, in units of cm)
  when \code{T = -2}.
}
\details{
Note that T (spectral period) can be a vector, while all other arguments
must be scalars.
}
\seealso{
  See \code{\link{Sa.nga}} for a comprehensive function that computes
  the ground motions from the AS08, BA08, CB08, and CY08 models, and
  outputs data in a matrix format.  See \code{\link{KBflatfile}} for an
  example of inputting and outputting earthquake data and predictions.
  
  For details on the sub procedures used for the individual NGA models, see
  \code{\link{subs.as}}, \code{\link{subs.ba}}, \code{\link{subs.cb}},
  and \code{\link{subs.cy}}.  See \code{\link{coefs}} for details on the
  period-independent model coefficients, and \code{\link{coefs.t.as}},
  \code{\link{coefs.t.ba}}, \code{\link{coefs.t.cb}}, and
  \code{\link{coefs.t.cy}} for details on the period-dependent model
  coefficients.
  
  For procedures on estimating input variables when they are not known,
  see \code{\link{Rx.calc}}, \code{\link{Rrup.calc}},
  \code{\link{dip.calc}}, \code{\link{W.calc}}, \code{\link{Ztor.calc}},
  \code{\link{Z1.calc}}, \code{\link{Z2.5.calc}}, and
  \code{\link{Zhyp.calc}}.
  These procedures are further described in Kaklamanos et al. (2011).

  For details on the spectral periods and ground motion parameters
  defined for each of the models, see \code{\link{modelPeriods}} or
  \code{\link{periods}}.  The functions \code{\link{getPeriod}} and
  \code{\link{interpolate}} provide methods for estimating spectral
  accelerations at intermediate periods between those with defined model
  coefficients.
}
\author{James Kaklamanos <james.kaklamanos@tufts.edu> and
  Eric M. Thompson <eric.thompson@tufts.edu>}
\references{
  
  Abrahamson, N., and W. Silva (2008). Summary of the Abrahamson & Silva
  NGA Ground-Motion Relations. \emph{Earthquake Spectra} \bold{24,}
  67--97.

  Atkinson, G. M., and D. M. Boore (2011). Modifications to Existing
  Ground-Motion Prediction Equations in Light of New
  Data. \emph{Bulletin of the Seismological Society of America},
  \bold{101,} 1121--1135.
  
  Boore, D. M., and G. M. Atkinson (2008). Ground-Motion Prediction
  Equations for the Average Horizontal Component of PGA, PGV, and
  5\%-Damped PSA at Spectral Periods between 0.01 s and 10.0
  s. \emph{Earthquake Spectra} \bold{24,} 99--138.

  Campbell, K. W., and Y. Bozorgnia (2007). Campbell-Bozorgnia NGA
  Ground Motion Relations for the Geometric Mean Horizontal Component of
  Peak and Spectral Ground Motion Parameters, \emph{PEER Report
    No. 2007/02}, Pacific Earthquake Engineering Research Center,
  University of California, Berkeley.
   
  Campbell, K. W., and Y. Bozorgnia (2008). NGA Ground Motion Model for
  the Geometric Mean Horizontal Component of PGA, PGV, PGD, and 5\%
  Damped Linear Elastic Response Spectra for Periods Ranging from 0.01
  to 10 s. \emph{Earthquake Spectra} \bold{24,} 139--171.

  Chiou, B. S.-J., and R. R. Youngs (2008). An NGA Model for the Average
  Horizontal Component of Peak Ground Motion and Response Spectra.
  \emph{Earthquake Spectra} \bold{24,} 173--215. 
  
  Kaklamanos, J., L. G. Baise, and D. M. Boore (2011).  Estimating
  Unknown Input Parameters when Implementing the NGA Ground-Motion
  Prediction Equations in Engineering Practice.
  \emph{Earthquake Spectra} \bold{27,} 1219--1235.
}
\examples{

# Assumed earthquake parameters for these examples:
M <- 7
Rjb <- 50
Rrup <- 51
Vs30 <- 300
Ztor <- 2
W <- 20
dip <- 80
VsFlag <- 0
Fhw <- 0
rake <- 180   # Strike-slip fault
Fas <- 0



####################################################################
# Example 1:  Illustration of the versatility of input for the Sa
#             functions (using CY08 as an example)

# Calculate PGA using the known input variables:
Sa.cy(M = M, Rjb = Rjb, Rrup = Rrup, Vs30 = Vs30, VsFlag = VsFlag,
      epsilon = 0, T = 0, dip = dip, W = W, Ztor = Ztor,
      rake = rake, Fhw = Fhw, AS = Fas)

# Alternately, the fault type may be input using the
# style-of-faulting flag variables:
Sa.cy(M = M, Rjb = Rjb, Rrup = Rrup, Vs30 = Vs30, VsFlag = VsFlag,
      epsilon = 0, T = 0, dip = dip, W = W, Ztor = Ztor,
      Frv = 0, Fnm = 0, Fhw = Fhw, AS = Fas)

# If the azimuth is known, it may be used in place of Fhw:
Sa.cy(M = M, Rjb = Rjb, Rrup = Rrup, Vs30 = Vs30, VsFlag = VsFlag,
      epsilon = 0, T = 0, dip = dip, W = W, Ztor = Ztor,
      Frv = 0, Fnm = 0, azimuth = -20, AS = Fas)

# The variables Rrup, dip, W, and Ztor may be left blank (or set
# to NA), and their defaults will be used in the calculation:
Sa.cy(M = M, Rjb = Rjb, Rrup = NA, Vs30 = Vs30, VsFlag = VsFlag,
      epsilon = 0, T = 0, Frv = 0, Fnm = 0, azimuth = -20, AS = Fas)



########################################################################
# Example 2:  Generate a plot of the predicted response spectrum (and
#             uncertainty) for a hypothetical earthquake using the AS08
#             model

# Redefine T to be a vector
# We only desire T >= 0.01 for plotting
T.list <- modelPeriods(model = "AS08", positive = TRUE) 

# Calculations
# Median
SaMedian <- Sa.as(M = M, Rjb = Rjb, Rrup = Rrup, Vs30 = Vs30,
                  VsFlag = VsFlag, T = T.list, dip = dip, W = W,
                  Ztor = Ztor, rake = rake, Fhw = Fhw, Fas = 0,
                  epsilon = 0)
# Median + 1 SD
SaPlus1SD <- Sa.as(M = M, Rjb = Rjb, Rrup = Rrup, Vs30 = Vs30,
                   VsFlag = VsFlag, T = T.list, dip = dip, W = W,
                   Ztor = Ztor, rake = rake, Fhw = Fhw, Fas = 0,
                   epsilon = 1)
# Median - 1 SD
SaMinus1SD <- Sa.as(M = M, Rjb = Rjb, Rrup = Rrup, Vs30 = Vs30,
                    VsFlag = VsFlag, T = T.list, dip = dip, W = W,
                    Ztor = Ztor, rake = rake, Fhw = Fhw, Fas = 0,
                    epsilon = -1)

# Plot
plot(T.list, SaMedian, type = "p", log = "xy", col = "blue", pch = 19,
     xlim = c(0.01, 10), ylim = c(0.001, 1), xaxs = "i", yaxs = "i",
     xlab = "Spectral Period, T [sec]", ylab = "Spectral Acceleration, Sa [g]",
     main = "AS08 Ground Motion Predictions:  Median +/- 1 SD")
points(T.list, SaMedian, pch = 19, col = "blue")
points(T.list, SaPlus1SD, pch = 19, col = "red")
points(T.list, SaMinus1SD, pch = 19, col = "red")
lines(T.list, SaMedian, lwd = 3, col = "blue")
lines(T.list, SaPlus1SD, lwd = 1, col = "red")
lines(T.list, SaMinus1SD, lwd = 1, col = "red")



#######################################################################
# Example 3:  Generate a plot of the median response spectra for the
#             same hypothetical earthquake, comparing the different
#             NGA models

# Calculations
# AS08
SaAS08 <- Sa.as(M = M, Rjb = Rjb, Rrup = Rrup, Vs30 = Vs30,
                VsFlag = VsFlag, epsilon = 0, T = T.list,
                dip = dip, W = W, Ztor = Ztor, rake = rake,
                Fhw = Fhw, Fas = 0)
# BA08
SaBA08 <- Sa.ba(M = M, Rjb = Rjb, Vs30 = Vs30, epsilon = 0,
                T = T.list, rake = rake)
# CB08
SaCB08 <- Sa.cb(M = M, Rjb = Rjb, Rrup = Rrup, Vs30 = Vs30,
                epsilon = 0, T = T.list, dip = dip, Ztor = Ztor,
                rake = rake)
# CY08
SaCY08 <- Sa.cy(M = M, Rjb = Rjb, Rrup = Rrup, Vs30 = Vs30,
                VsFlag = VsFlag, epsilon = 0, T = T.list,
                dip = dip, W = W, Ztor = Ztor, rake = rake,
                Fhw = Fhw, AS = 0)

# Plot
plot(T.list, SaAS08, type = "l", log = "xy", col = "blue", pch = 19, lwd = 2,
     xlim = c(0.01, 10), ylim = c(0.001, 1), xaxs = "i", yaxs = "i",
     xlab = "Spectral Period, T [sec]", ylab = "Spectral Acceleration, Sa [g]",
     main = "Comparison of NGA Ground Motion Predictions")
lines(T.list, SaBA08, lwd = 2, col = "red")
lines(T.list, SaCB08, lwd = 2, col = "darkgreen")
lines(T.list, SaCY08, lwd = 2, col = "black")
legend(x = "bottomleft", inset = 0.02, lwd = 2, bty = "n",
       col = c("blue", "red", "darkgreen", "black"),
       legend = c("AS08", "BA08", "CB08", "CY08"))
}
