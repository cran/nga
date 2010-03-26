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
  Generation Attenuation of Ground Motions (NGA) project in 2008.
}
\usage{
Sa.as(M, Rjb, Rrup = NA, Rx = NA, rake = NA, Frv = NA, Fnm = NA,
      Fhw = NA, dip = NA, W = NA, Ztor = NA, Vs30, Z1.0 = NA,
      VsFlag, Fas = 0, Zhyp = NA, azimuth = NA, epsilon, T)
Sa.ba(M, Rjb, Vs30, rake = NA, U = NA, SS = NA, NS = NA, RS = NA,
      epsilon, T)
Sa.cb(M, Rjb, Rrup = NA, rake = NA, Frv = NA, Fnm = NA, dip = NA,
      W = NA, Ztor = NA, Vs30, Z1.0 = NA, Z1.5 = NA, Z2.5 = NA,
      Zhyp = NA, Fhw = NA, azimuth = NA, arb = 0, epsilon, T)
Sa.cy(M, Rjb, Rrup = NA, Rx = NA, rake = NA, Frv = NA, Fnm = NA,
      Fhw = NA, dip = NA, W = NA, Ztor = NA, Vs30, Z1.0 = NA,
      VsFlag, AS = 0, Zhyp = NA, azimuth = NA, epsilon, T)
}
\arguments{
  \item{M}{Moment magnitude of earthquake.}
  \item{Rjb}{Joyner-Boore distance (km):  the horizontal distance to the
    surface projection of the rupture plane.}
  \item{Rrup}{Rupture distance (km):  the closest distance to the
    rupture plane; if left empty (or set negative), Rrup is calculated
    from Rx, the source-to-site azimuth, and the geometric rupture
    parameters (Ztor, W, and dip) using \code{\link{Rrup.calc}}.}
  \item{Rx}{Site coordinate (km):  The horizontal distance to the
    surface projection of the top edge of the rupture plane, measured
    perpendicular to the strike.  If left empty, Rx is calculated from Rjb, the
    source-to-site azimuth, and the geometric rupture parameters (Ztor,
    W, and dip) using \code{\link{Rx.calc}}.  When only Rjb and the
    azimuth are assumed, Rjb is used to calculate Rx, which is then used
    to calculate Rrup.}
  \item{rake}{Rake angle of fault movement (deg).  Either the rake angle
    or the style-of-faulting flag variables (Frv and Fnm for AS08, CB08,
    and CY08; and U, SS, NS, and SS for BA08) must be specified.}
  \item{Frv}{Reverse style-of-faulting flag (1 for reverse faulting, 0
    otherwise) for the AS08, CB08, and CY08 models.  Either (a) the rake
    angle, or (b) both Frv and Fnm, must be specified.  Reverse faulting
    is characterized by rake angles in the range 30 <= rake <= 150 deg.}
  \item{Fnm}{Normal style-of-faulting flag (1 for normal faulting, 0
    otherwise) for the AS08, CB08, and CY08 models.  Either (a) the rake
    angle, or (b) both Frv and Fnm, must be specified.  Normal faulting
    is characterized by rake angles in the range -120 <= rake <= -60 deg for
    the AS08 and CY08 models and -150 <= rake <= -30 deg for the CB08
    model.}
  \item{U}{Unspecified style-of-faulting flag for the BA08 model;
    equal to 1 if the user wishes to perform a generic ground motion
    calculation when the style of faulting is unspecified, and 0 otherwise.} 
  \item{RS}{Reverse style-of-faulting flag for the BA08 model;
    equal to 1 for reverse faulting (30 <= rake <= 150 deg), and 0 otherwise.}
  \item{NS}{Normal style-of-faulting flag for the BA08 model;
    equal to 1 for normal faulting (-150 <= rake <= -30 deg), and 0 otherwise.}
  \item{SS}{Strike-slip style-of-faulting flag for the BA08 model; equal to 1 for
    strike-slip faulting (when the rake is not in either of the ranges
    specified for RS or NS), and 0 otherwise.}
  \item{dip}{Dip angle of the rupture plane (deg).  If left empty (or set
    negative), the dip is estimated using \code{\link{dip.calc}}.}
  \item{W}{Down-dip width of rupture plane (km).  If left empty (or set
    negative), W is estimated using \code{\link{W.calc}}.}
  \item{Ztor}{Depth to top of rupture (km).  If left empty (or set
    negative), Ztor is estimated using \code{\link{Ztor.calc}}.}
  \item{Vs30}{Time-averaged shear wave velocity over a subsurface
    depth of 30 meters (m/s).}
  \item{VsFlag}{Flag variable indicating how Vs30 is obtained; equal
    to 1 if Vs30 is measured, and 0 if Vs30 is estimated or inferred.}
  \item{Z1.0}{Depth to Vs = 1.0 km/s (m).  If left empty (or set
    negative), Z1.0 is estimated using \code{\link{Z1.calc.as}} for the
    AS08 model and \code{\link{Z1.calc.cy}} for the CY08 model.}
  \item{Z1.5}{Depth to Vs = 1.5 km/s (m).  Z1.5 is not utilized in
    ground motion calculations, but if available, it may be used to
    estimate Z2.5 for the CB08 model.}
  \item{Z2.5}{Depth to Vs = 2.5 km/s (km).  If left empty (or set negative),
    Z2.5 is estimated from Z1.5 or Z1.0 if available, using the
    recommendations in Campbell and Bozorgnia (2007).  If neither Z1.5 nor
    Z1.0 is available, then Vs30 is used to estimate Z1.0 using
    \code{\link{Z1.calc.as}}, which is in turn used to estimate Z2.5.}
  \item{Fas}{Aftershock flag for AS08; equal to 1 for aftershocks and 0 for
    mainshocks (the default).}
  \item{AS}{Aftershock flag for CY08; equal to 1 for aftershocks and 0 for
    mainshocks (the default).}
  \item{Fhw}{Hanging wall flag; equal to 1 for sites on the hanging
    wall side of the fault (Rx >= 0; azimuth >= 0), and 0 otherwise.}
  \item{Zhyp}{Hypocentral depth of the earthquake (km).  Zhyp is not
    utilized in ground motion calculations, but it may be used to estimate
    Ztor.  See \code{\link{Ztor.calc}} for details.}
  \item{azimuth}{Source-to-site azimuth (deg); see Figure 2 in
    Kaklamanos and Baise (2010) for details.  Used by
    \code{\link{Rx.calc}} and \code{\link{Rrup.calc}} for distance
    calculations.}
  \item{arb}{Flag variable indicating the method of determining aleatory
    uncertainty for the CB08 model; equal to 1 if the standard deviation
    should be calculated for the arbitrary horizontal component of ground
    motion, and 0 if the standard deviation should be calculated for the
    geometric mean horizontal ground motion (the default).}
  \item{epsilon}{number of standard deviations to be considered in the
    calculations.  Use 0 to obtain a median estimate of ground motion.}
  \item{T}{Spectral period (sec).  Use 0 for PGA, -1 for PGV, and -2 for
    PGD (CB08 only).  For spectral acceleration, T must be in the range
    0.01 <= T <= 10 sec.  If the specified period is within the
    allowable range and does not have defined equations, the program
    uses log-log interpolation (using \code{\link{interpolate}}) between
    the next-highest and next-lowest spectral periods with defined
    equations.}
}
\value{
  The spectral acceleration (in units of \emph{g}) at period T; peak
  ground acceleration (PGA, in units of \emph{g}) when T = 0; peak
  ground velocity (PGV, in units of cm/sec) when T = -1; and peak ground
  displacement using the CB08 model (PGD, in units of cm) when T = -2
  sec.
}
\details{
Note that T (spectral period) can be a vector, while all other arguments are 
scalars.
}
\seealso{
  See \code{\link{Sa.nga}} for a comprehensive function that computes the ground
  motions from the AS08, BA08, CB08, and CY08 models, and outputs data
  in a matrix format.  See \code{\link{mainshocks}} for an example of
  inputting and outputting earthquake data and predictions.
  
  For details on the sub procedures used for the individual NGA models, see
  \code{\link{subs.as}}, \code{\link{subs.ba}}, \code{\link{subs.cb}},
  and \code{\link{subs.cy}}.  See \code{\link{coefs}} for details on the
  period-independent model coefficients, and \code{\link{coefs.t.as}},
  \code{\link{coefs.t.ba}}, \code{\link{coefs.t.cb}}, and
  \code{\link{coefs.t.cy}} for details on the period-dependent model coefficients.
  
  For procedures on estimating input variables when they are not known,
  see \code{\link{Rx.calc}}, \code{\link{Rrup.calc}},
  \code{\link{dip.calc}}, \code{\link{W.calc}}, \code{\link{Ztor.calc}},
  and \code{\link{Z1.calc}}.

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

  Kaklamanos, J., and L. G. Baise (2010).  Relationships between
  Distance Measures in Recent Ground Motion Prediction Equations.
  \emph{Seismological Research Letters} (in review).
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



# Example 1:  Illustration of the versatility of input for the Sa
#             functions (using CY08 as an example)

# Calculate PGA using the known input variables:
Sa.cy(M = M, Rjb = Rjb, Rrup = Rrup, rake = rake, Fhw = Fhw,
      dip = dip,  W = W, Ztor = Ztor, Vs30 = Vs30,
      VsFlag = VsFlag, AS = Fas, epsilon = 0, T = 0)

# Alternately, the fault type may be input using the
# style-of-faulting flag variables:
Sa.cy(M = M, Rjb = Rjb, Rrup = Rrup, Frv = 0, Fnm = 0, Fhw = Fhw,
      dip = dip,  W = W, Ztor = Ztor, Vs30 = Vs30,
      VsFlag = VsFlag, AS = Fas, epsilon = 0, T = 0)

# If the azimuth is known, it may be used in place of Fhw:
Sa.cy(M = M, Rjb = Rjb, Rrup = Rrup, Frv = 0, Fnm = 0, azimuth = -20,
      dip = dip,  W = W, Ztor = Ztor, Vs30 = Vs30,
      VsFlag = VsFlag, AS = Fas, epsilon = 0, T = 0)

# The variables Rrup, dip, W, and Ztor may be left blank, and their
# defaults will be used in the calculation:
Sa.cy(M = M, Rjb = Rjb, Frv = 0, Fnm = 0, azimuth = -20,
      Vs30 = Vs30, VsFlag = VsFlag, AS = Fas, epsilon = 0, T = 0)



# Example 2:  Generate a plot of the predicted response spectrum (and
#             uncertainty) for a hypothetical earthquake using the AS08
#             model

# Redefine T to be a vector
# We only desire T >= 0.01 for plotting
T.list <- modelPeriods(model = "AS08", positive = TRUE) 

# Calculations
# Median
SaMedian <- Sa.as(M = M, Rjb = Rjb, Rrup = Rrup, rake = rake,
                  Fhw = Fhw, dip = dip, W = W, Ztor = Ztor, Vs30 = Vs30,
                  VsFlag = VsFlag, Fas = 0, epsilon = 0, T = T.list)
# Median + 1 SD
SaPlus1SD <- Sa.as(M = M, Rjb = Rjb, Rrup = Rrup, rake = rake,
                   Fhw = Fhw, dip = dip, W = W, Ztor = Ztor, Vs30 = Vs30,
                   VsFlag = VsFlag, Fas = 0, epsilon = 1, T = T.list)
# Median - 1 SD
SaMinus1SD <- Sa.as(M = M, Rjb = Rjb, Rrup = Rrup, rake = rake,
                    Fhw = Fhw, dip = dip, W = W, Ztor = Ztor, Vs30 = Vs30,
                    VsFlag = VsFlag, Fas = 0, epsilon = -1, T = T.list)

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



# Example 3:  Generate a plot of the median response spectra for the
#             same hypothetical earthquake, comparing the different
#             NGA models

# Calculations
# AS08
SaAS08 <- Sa.as(M = M, Rjb = Rjb, Rrup = Rrup, rake = rake,
                Fhw = Fhw, dip = dip, W = W, Ztor = Ztor, Vs30 = Vs30,
                VsFlag = VsFlag, Fas = 0, epsilon = 0, T = T.list)
# BA08
SaBA08 <- Sa.ba(M = M, Rjb = Rjb, Vs30 = Vs30, rake = rake,
                epsilon = 0, T = T.list)
# CB08
SaCB08 <- Sa.cb(M = M, Rjb = Rjb, Rrup = Rrup, rake = rake, dip = dip,
                Ztor = Ztor, Vs30 = Vs30, epsilon = 0, T = T.list)
# CY08
SaCY08 <- Sa.cy(M = M, Rjb = Rjb, Rrup = Rrup, rake = rake, Fhw = Fhw,
                dip = dip,  W = W, Ztor = Ztor, Vs30 = Vs30,
                VsFlag = VsFlag, AS = Fas, epsilon = 0, T = T.list)

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
