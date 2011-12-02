\name{Ground Motion Predictions for all NGA Models}
\alias{Sa.nga}
\alias{nga}
\title{Ground Motion Predictions for all NGA Models}
\description{
  Comprehensive function that estimates ground motion parameters using the
  AS08, BA08, CB08, and CY08 models from the Next Generation Attenuation
  of Ground Motions (NGA) project in 2008.  The function \code{Sa.nga}
  is designed to mimic the output from Boore and Campbell's Fortran
  output files.
}
\usage{
Sa.nga(M, Rjb, Vs30, T, Rrup = NA, Rx = NA, dip = NA, W = NA,
       Ztor = NA, Z1.0 = NA, Z1.5 = NA, Z2.5 = NA, rake = NA,
       Frv = NA, Fnm = NA, Fhw = NA, azimuth = NA, Zhyp = NA,
       Fas = 0, epsilon = 1)
}
\arguments{
  \item{M}{Moment magnitude of earthquake.}
  \item{Rjb}{Joyner-Boore distance (km):  the horizontal distance to the
    surface projection of the rupture plane.}
  \item{Vs30}{Time-averaged shear wave velocity over a subsurface
    depth of 30 meters (m/s).}
  \item{T}{Spectral period (sec).  Use 0 for PGA and -1 for PGV.  For
    spectral acceleration, \code{T} must be in the range \code{0.01 <= T
      <= 10 sec}.  If the  specified period is within the allowable
    range and does not have defined equations, the program uses log-log
    interpolation (using \code{\link{interpolate}}) between the
    next-highest and next-lowest spectral periods with defined equations.} 
  \item{Rrup}{Rupture distance (km):  the closest distance to the
    rupture plane; if left empty, \code{Rrup} is calculated from
    \code{Rx}, the source-to-site azimuth, and the geometric rupture
    parameters (\code{Ztor}, \code{W}, and \code{dip}) using
    \code{\link{Rrup.calc}}.}
  \item{Rx}{Site coordinate (km):  The horizontal distance to the
    surface projection of the top edge of the rupture plane, measured
    perpendicular to the strike.  If left empty, \code{Rx} is calculated from
    \code{Rjb}, the source-to-site azimuth, and the geometric rupture
    parameters (\code{Ztor}, \code{W}, and \code{dip}) using
    \code{\link{Rx.calc}}.  When only \code{Rjb} and the azimuth are
    assumed, \code{Rjb} is used to calculate \code{Rx},
    which is then used to calculate \code{Rrup}.}
  \item{dip}{Dip angle of the rupture plane (deg).  If left empty, the
    dip is estimated using \code{\link{dip.calc}}.}
  \item{W}{Down-dip width of rupture plane (km).  If left empty, \code{W} is
    estimated using \code{\link{W.calc}}.}
  \item{Ztor}{Depth to top of rupture (km).  If left empty, \code{Ztor} is
    estimated using \code{\link{Ztor.calc}}.}
  \item{Z1.0}{Depth to Vs = 1.0 km/s (m).  If left empty, \code{Z1.0} is
    estimated using \code{\link{Z1.calc.as}} for the AS08 model and
    \code{\link{Z1.calc.cy}} for the CY08 model.}
  \item{Z1.5}{Depth to Vs = 1.5 km/s (m).  \code{Z1.5} is not utilized in
    ground motion calculations, but if available, it may be used to
    estimate \code{Z2.5} for the CB08 model.}
  \item{Z2.5}{Depth to Vs = 2.5 km/s (m; \bold{note the units}).  If
    left empty, \code{Z2.5} is estimated from \code{Z1.5} or \code{Z1.0}
    if available, using the recommendations in Campbell and Bozorgnia
    (2007); see \code{\link{Z2.5.calc}}.  If neither \code{Z1.5} nor
    \code{Z1.0} is available, then \code{Vs30} is used to estimate
    \code{Z1.0} using \code{\link{Z1.calc.as}}, which is in turn used to
    estimate \code{Z2.5}.}
  \item{rake}{Rake angle of fault movement (deg).  Either the rake angle
    or the style-of-faulting flag variables (\code{Frv} and \code{Fnm})
    must be specified.}
  \item{Frv}{Reverse style-of-faulting flag (1 for reverse faulting, 0
    otherwise).  Either (a) the rake angle, or (b) both \code{Frv} and
    \code{Fnm}, must be specified.  Reverse faulting is characterized by
    rake angles in the range \code{30 <= rake <= 150 deg} for the AS08,
    BA08, and CY08 models; and in the range \code{30 < rake < 150 deg}
    for the CB08 model.}
  \item{Fnm}{Normal style-of-faulting flag (1 for normal faulting, 0
    otherwise).  Either (a) the rake angle, or (b) both \code{Frv} and
    \code{Fnm}, must be specified.  Normal faulting is characterized by
    rake angles in the range \code{-120 <= rake <= -60 deg} for the AS08
    and CY08 models, \code{-150 <= rake <= -30 deg} for the BA08 model,
    and \code{-150 < rake < -30 deg} for the CB08 model.}
  \item{Fhw}{Hanging wall flag; equal to 1 for sites on the hanging
    wall side of the fault (\code{Rx >= 0; azimuth >= 0)}, and 0
    otherwise.  Either \code{Fhw}, \code{Rx}, or the azimuth must be
    specified.}
  \item{azimuth}{Source-to-site azimuth (deg); see Kaklamanos et al.
    (2011) for details.  Used by
    \code{\link{Rx.calc}} and \code{\link{Rrup.calc}} for distance
    calculations.  Either \code{Fhw}, \code{Rx}, or the azimuth must be
    specified.}
  \item{Zhyp}{Hypocentral depth of the earthquake (km).  \code{Zhyp} is not
    utilized in ground motion calculations, but it may be used to estimate
    \code{Ztor}.  See \code{\link{Ztor.calc}} for details.}
  \item{Fas}{Aftershock flag; equal to 1 for aftershocks and 0 for
    mainshocks (the default)}
  \item{epsilon}{Number of standard deviations to be considered in the
    calculations (default value is 1). The function \code{Sa.nga}
    automatically outputs the median estimates (corresponding to
    \code{epsilon = 0}) as well as the estimates corresponding to the
    median estimate plus and minus \code{epsilon * sigmaTotal}}
    
}
\value{
  The function \code{Sa.nga} outputs a data frame composed of the
  following 62 columns:

  \bold{Input Variables:}
  \tabular{ll}{
    \code{T} \tab Spectral period, sec [input] \cr
    \code{M} \tab Moment magnitude [input] \cr
    \code{Rjb} \tab Joyner-Boore distance (km) [input] \cr
    \code{Rrup.in} \tab Rupture distance (km) [input] \cr
    \code{Rrup.out} \tab Rupture distance (km) [calculated if
    \code{Rrup.in} is not specified] \cr
    \code{Rx.in} \tab Site coordinate (km) [input] \cr
    \code{Rx.out} \tab Site coordinate (km) [calculated if \code{Rx.in} is not
    specified] \cr
    \code{azimuth.in} \tab source-to-site azimuth (deg) [input] \cr
    \code{azimuth.out} \tab source-to-site azimuth (deg) [calculated if
    \code{azimuth.in} is not specified] \cr
    \code{Fhw} \tab Hanging wall flag \cr
    \code{Zhyp.in} \tab hypocentral depth (km) [input] \cr
    \code{Zhyp.out} \tab hypocentral depth (km) [calculated if
    \code{Zhyp.in} is not specified] \cr
    \code{rake.in} \tab Rake angle of fault movement (deg) [input] \cr
    \code{rake.out} \tab Rake angle of fault movement (deg) [calculated
    if \code{rake.in} is not specified] \cr
    \code{Frv1} \tab Reverse style-of-faulting flag for AS08, BA08, and
    CY08 [input] \cr
    \code{Frv2.cb} \tab Reverse style-of-faulting flag for CB08 \cr
    \code{Fnm1} \tab Normal style-of-faulting flag for AB08 and CY08 \cr
    \code{Fnm2.ba} \tab Normal style-of-faulting flag for BA08 \cr
    \code{Fnm3.cb} \tab Normal style-of-faulting flag for CB08 \cr
    \code{dip.in} \tab Fault dip angle (deg) [input] \cr
    \code{dip.out} \tab Fault dip angle (deg) [calculated if \code{dip.in} is
    not specified] \cr
    \code{W.in} \tab Down-dip rupture width (km) [input] \cr
    \code{W.out} \tab Down-dip rupture width (km) [calculated if
    \code{W.in} is not specified] \cr
    \code{Ztor.in} \tab Depth to top of rupture (km) [input] \cr
    \code{Ztor.out} \tab Depth to top of rupture (km) [calculated if
    \code{Ztor.in} is not specified] \cr
    \code{Vs30} \tab Time-averaged shear wave velocity over 30 m
    subsurface depth  (m/sec) [input] \cr
    \code{Z1.0in} \tab Depth to Vs of 1.0 km/sec  (m) [input] \cr
    \code{Z1.0as} \tab Depth to Vs of 1.0 km/sec  (m) [calculated for
    use in AS08 model] \cr
    \code{Z1.0cy} \tab Depth to Vs of 1.0 km/sec  (m) [calculated for
    use in CY08 model] \cr
    \code{Z1.5in} \tab Depth to Vs of 1.5 km/sec  (m) [input] \cr
    \code{Z2.5in} \tab Depth to Vs of 2.5 km/sec  (m) [input] \cr
    \code{Z2.5out} \tab Depth to Vs of 2.5 km/sec  (m) [calculated from
    \code{Z1.0} for use in CB08 model] \cr
    \code{Fas} \tab Aftershock flag [input] \cr
    \code{epsilon} \tab number of standard deviations considered in the
    calculations [input] \cr
  }
  \bold{Output Variables:}
  \tabular{ll}{
    \bold{AS08 Model:} \tab  \cr
    \code{Y50.as} \tab Median ground motion estimate using AS08 (\code{epsilon}
    = 0) \cr
    \code{YplusEpsilon.meas.as} \tab Upper ground motion estimate
    using AS08, for measured \code{Vs30} (\code{VsFlag = 1}) \cr
    \code{YplusEpsilon.est.as} \tab Upper ground motion estimate
    using AS08, for estimated \code{Vs30} (\code{VsFlag = 0}) \cr
    \code{YminusEpsilon.meas.as} \tab Lower ground motion estimate
    using AS08, for measured \code{Vs30} (\code{VsFlag = 1}) \cr
    \code{YminusEpsilon.est.as} \tab Lower ground motion estimate
    using AS08, for estimated \code{Vs30} (\code{VsFlag = 0}) \cr
    \code{sdMeas.as} \tab total standard deviation using AS08, for
    measured \code{Vs30} (\code{VsFlag = 1}) \cr
    \code{sdEst.as} \tab total standard deviation using AS08, for
    estimated \code{Vs30} (\code{VsFlag = 0}) \cr
    \bold{BA08 Model:} \tab  \cr
    \code{Y50.ba} \tab Median ground motion estimate using BA08 \cr
    \code{Y50mod.ba} \tab Median ground motion estimate using modified
    BA08 (\code{AB11 = 1}) \cr
    \code{YplusEpsilon.ba} \tab Upper ground motion estimate
    using BA08 \cr
    \code{YplusEpsilon.mod.ba} \tab Upper ground motion estimate
    using modified BA08 (\code{AB11 = 1}) \cr
    \code{YminusEpsilon.ba} \tab Lower ground motion estimate
    using BA08 \cr
    \code{YminusEpsilon.mod.ba} \tab Lower ground motion estimate
    using modified BA08 (\code{AB11 = 1}) \cr
    \code{sd.ba} \tab total standard deviation using BA08 \cr
    \bold{CB08 Model:} \tab  \cr
    \code{Y50.cb} \tab Median ground motion estimate using CB08 (\code{epsilon}
    = 0) \cr
    \code{YplusEpsilon.GM.cb} \tab Upper CB08 estimate for the
    geometric mean horizontal component (\code{arb = 0}) \cr
    \code{YplusEpsilon.arb.cb} \tab Upper CB08 estimate for the
    arbitrary horizontal component (\code{arb = 1}) \cr
    \code{YminusEpsilon.GM.cb} \tab Lower CB08 estimate for the
    geometric mean horizontal component (\code{arb = 0}) \cr
    \code{YminusEpsilon.arb.cb} \tab Lower CB08 estimate for the
    arbitrary horizontal component (\code{arb = 1}) \cr
    \code{sdGM.cb} \tab CB08 total standard deviation
    for the geometric mean horizontal component (\code{arb = 0}) \cr
    \code{sdArb.cb} \tab CB08 total standard deviation for the
    arbitrary horizontal component (\code{arb = 1}) \cr
    \bold{CY08 Model:} \tab  \cr
    \code{Y50.cy} \tab Median ground motion estimate using CY08 (\code{epsilon}
    = 0) \cr
    \code{YplusEpsilon.meas.cy} \tab Upper ground motion estimate
    using CY08, for measured \code{Vs30} (\code{VsFlag = 1}) \cr
    \code{YplusEpsilon.est.cy} \tab Upper ground motion estimate
    using CY08, for estimated \code{Vs30} (\code{VsFlag = 0}) \cr
    \code{YminusEpsilon.meas.cy} \tab Lower ground motion estimate
    using CY08, for measured \code{Vs30} (\code{VsFlag = 1}) \cr
    \code{YminusEpsilon.est.cy} \tab Lower ground motion estimate
    using CY08, for estimated \code{Vs30} (\code{VsFlag = 0}) \cr
    \code{sdMeas.cy} \tab total standard deviation using CY08, for
    measured \code{Vs30} (\code{VsFlag = 1}) \cr
    \code{sdEst.cy} \tab total standard deviation using CY08, for
    estimated \code{Vs30} (\code{VsFlag = 0}) \cr
     \tab
  }

}
\details{
Note that T (spectral period) can be a vector, while all other arguments are 
scalars.

In the \dQuote{Output Section} of this function, \dQuote{Y} refers to
the ground motion parameter of interest, which can be:
\tabular{l}{
  Sa = Spectral acceleration (g); \cr
  PGA = Peak ground acceleration (g), calculated by evaluating Sa at
  \code{T = 0}; \cr
  PGV = Peak ground velocity (cm/sec), calculated by evaluating Sa at
  \code{T = -1}.
}
Because only the CB08 model has coefficients for PGD (peak ground
displacement), the CB08-specific function \code{\link{Sa.cb}} must be
used to obtain predictions for PGD.  In addition, \dQuote{sd} refers to
the standard deviation of the ground motion estimate, which is presented
in natural log space.

The flag variables \code{VsFlag} and \code{arb} refer to:

\code{VsFlag} = Flag variable indicating how \code{Vs30} is obtained
(AS08 and CY08 models only); equal to 1 if \code{Vs30} is measured,
and 0 if \code{Vs30} is estimated or inferred.

\code{arb} = Flag variable indicating the method of determining aleatory
uncertainty for the CB08 model; equal to 1 if the standard deviation
should be calculated for the arbitrary horizontal component of ground
motion, and 0 if the standard deviation should be calculated for the
geometric mean horizontal ground motion.

These two indicator variables represent model-specific options for
output: AS08 and CY08 have different standard deviation terms for
measured and inferred Vs30 (specified by \code{VsFlag}), and CB08 is the
only model that offers predictions for the arbitrary horizontal
component of ground motion (\code{arb}).  For each case (0 and 1) of
each of these three indicator variables, \code{Sa.nga} provides the
estimated ground motion.  This output is consistent with that of the
Fortran program described later in this report.
The model-specific functions \code{\link{Sa.as}}, \code{\link{Sa.ba}},
\code{\link{Sa.cb}}, and \code{\link{Sa.cy}} allow the user to specify
the values of the indicator variables in the arguments to the functions.

The median BA08 estimate is presented in terms of the original GMPE
(Boore and Atkinson, 2008) as well as the modified GMPE given by
Atkinson and Boore (2011).  The small-magnitude modification affects
ground motion estimates for \code{M <= 5.75}.  The modified BA08 model
corresponds to \code{AB11 = 1} in the \code{\link{Sa.ba}} function, and
the original BA08 model corresponds to \code{AB11 = 0}.

}
\seealso{
  See \code{\link{Sa.as}}, \code{\link{Sa.ba}}, \code{\link{Sa.cb}}, and
  \code{\link{Sa.cy}} for separate functions that compute ground motion
  parameters using the individual NGA models. See
  \code{\link{KBflatfile}} for an example of inputting and outputting
  earthquake data and predictions.
    
  For details on the sub procedures used for the individual NGA models,
  see \code{\link{subs.as}}, \code{\link{subs.ba}},
  \code{\link{subs.cb}}, and \code{\link{subs.cy}}.  See
  \code{\link{coefs}} for details on the period-independent model
  coefficients, and \code{\link{coefs.t.as}}, \code{\link{coefs.t.ba}},
  \code{\link{coefs.t.cb}}, and \code{\link{coefs.t.cy}} for details on
  the period-dependent model coefficients.
   
  For procedures on estimating input variables when they are not known,
  see \code{\link{Rx.calc}}, \code{\link{Rrup.calc}},
  \code{\link{dip.calc}}, \code{\link{W.calc}}, \code{\link{Ztor.calc}},
  \code{\link{Z1.calc}}, \code{\link{Z2.5.calc}}, and
  \code{\link{Zhyp.calc}}.  These procedures are further described in
  Kaklamanos et al. (2011).

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
  Data. \emph{Bulletin of the Seismological Society of America}
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
T.list <- c(0, 0.1, 0.5, 1)
dip <- 80
W <- 20
Ztor <- 2
rake <- 180
Fhw <- 0
Fas <- 0




#########################################################
# Example 1:  Illustration of the versatility of input

# First calculate ground motions using the known input variables
# Some of the variables (such as Z1.0 and Rx) are unknown, and will
# be calculated by the program
Sa.nga(M = M, Rjb = Rjb, Vs30 = Vs30, epsilon = 1, T = T.list,
       Rrup = Rrup, dip = dip, W = W, Ztor = Ztor, rake = rake,
       Fhw = Fhw, Fas = Fas)

# Repeat the ground motion calculation the bare minimum necessary requirements
Sa.nga(M = M, Rjb = Rjb, Vs30 = Vs30, epsilon = 1, T = T.list,
       rake = rake, Fhw = Fhw)

# Note that the style-of-faulting flag variables may be used in place
# of the rake, and that the azimuth (if known) may be used instead of Fhw
Sa.nga(M = M, Rjb = Rjb, Vs30 = Vs30, epsilon = 1, T = T.list,
       Frv = 0, Fnm = 0, azimuth = -30)



#######################################################################
# Example 2: Generate a plot of the predicted response spectrum (and
#            uncertainty) for a hypothetical earthquake using the BA08
#            model

# Redefine T to be a vector
# We only desire T >= 0.01 for plotting
T.list <- modelPeriods(model = "BA08", positive = TRUE) 

# Ground motion calculations
ResultsMatrix <- Sa.nga(M = M, Rjb = Rjb, Rrup = Rrup, Vs30 = Vs30,
                        epsilon = 1, T = T.list, dip = dip, W = W,
                        Ztor = Ztor, rake = rake, Fhw = Fhw, Fas = Fas)

# To see the names of all the columns in the data frame, use the "names"
# function on a column of the matrix:
names(ResultsMatrix)

# To access individual columns of the data frame, use the "$" operator:
SaMedian <- ResultsMatrix$Y50.ba
SaPlusEpsilon <- ResultsMatrix$YplusEpsilon.ba
SaMinusEpsilon <- ResultsMatrix$YminusEpsilon.ba

# Plot
plot(T.list, SaMedian, type = "p", log = "xy", col = "blue", pch = 19,
     xlim = c(0.01, 10), ylim = c(0.001, 1), xaxs = "i", yaxs = "i",
     xlab = "Spectral Period, T [sec]", ylab = "Spectral Acceleration, Sa [g]",
     main = "BA08 Ground Motion Predictions:  Median +/- 1 SD")
points(T.list, SaMedian, pch = 19, col = "blue")
points(T.list, SaPlusEpsilon, pch = 19, col = "red")
points(T.list, SaMinusEpsilon, pch = 19, col = "red")
lines(T.list, SaMedian, lwd = 3, col = "blue")
lines(T.list, SaPlusEpsilon, lwd = 1, col = "red")
lines(T.list, SaMinusEpsilon, lwd = 1, col = "red")



#######################################################################
# Example 3:  Generate a plot of the median response spectra for the
#             same hypothetical earthquake, comparing the different
#             NGA models

# Redefine T to be a vector
# We only desire T >= 0.01 for plotting
T.list <- modelPeriods(model = "BA08", positive = TRUE) 

# Ground motion calculations
ResultsMatrix <- Sa.nga(M = M, Rjb = Rjb, Rrup = Rrup, Vs30 = Vs30,
                        epsilon = 1, T = T.list, dip = dip, W = W,
                        Ztor = Ztor, rake = rake, Fhw = Fhw, Fas = Fas)

# Access individual columns of the data frame using the "$" operator:
SaAS08 <- ResultsMatrix$Y50.as
SaBA08 <- ResultsMatrix$Y50.ba
SaCB08 <- ResultsMatrix$Y50.cb
SaCY08 <- ResultsMatrix$Y50.cy

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
