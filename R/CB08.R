# CAMPBELL & BOZORGNIA NGA MODEL
# Campbell, K. W., and Y. Bozorgnia. NGA Ground Motion Model for Geometric Mean
# Horizontal Component of PGA, PGV, PGD, and 5% Damped Linear Elastic Response
# Spectra for Periods Ranging from 0.01 to 10 s. Earthquake Spectra, Vol. 24, pp. 139-171.


# OUTLINE OF CODE
# 1. Model Coefficients
#    a. Periods with defined coefficients
#    b. Period-independent coefficients
#    c. Period-dependent coefficients for median ground motion term
#    d. Coefficients for standard deviation term
# 2. Necessary Functions for Calculating Median Ground Motion
#    a. Magnitude Term
#    b. Distance Term
#    c. Fault Mechanism Term
#    d. Hanging-Wall Term
#    e. Shallow Site Response Term
#    f. Basin Response Term
#    g. Calculation of A1100
# 3. Necessary Functions for Calculating Standard Deviation Term
#    a. Partial derivative of Fsite with respect to A1100
#    b. Sigma, intra-event standard deviation
#    c. Tau, inter-event standard deviation
#    d. SigmaTot, total standard deviation (geometric mean)
#    e. SigmaArb, total standard deviation (arbitrary horizontal component)
# 4. Median Ground Motion Calculation of Sa, PGA, and PGV
# 5. Final Function for CB08 Ground Motion Calculations
#    a. Check input parameters
#    b. Obtain estimates of unspecified input parameters
#    c. Calculate ground motion parameter
#    **** This is the function that users will primarily use for ground motion calculations ****


# INPUT PARAMETERS FOR THE FUNCTIONS BELOW:
#   M = Moment magnitude
#   Rjb = Joyner-Boore distance (km)
#   Rrup = Rupture distance (km)
#   Rx = Site coordinate (km); used for the calculation of Rrup
#   rake = Rake angle of fault movement (deg)
#   dip = Fault dip (deg)
#   W = Down-dip rupture width (km); used for the calculation of Rrup
#   Ztor = Depth to top of rupture (km)
#   Vs30 = Time-averaged shear wave velocity over 30 m subsurface depth  (m/sec)
#   Z1.0 = Depth to Vs = 1.0 km/sec  (m); used to estimate Z2.5
#   Z1.5 = Depth to Vs = 1.5 km/sec  (m); used to estimate Z2.5
#   Z2.5 = Depth to Vs = 2.5 km/sec  (km) -- note the units difference from Z1.0 and Z1.5
#   Frv = Reverse style-of-faulting flag (1 for reverse faulting,
#         0 otherwise); calculated from rake
#   Fnm = Normal style-of-faulting flag (1 for normal faulting,
#         0 otherwise); calculated from rake
#   Zhyp = hypocentral depth (km)
#   azimuth = source-to-site azimuth (deg); see Figure 2 in Kaklamanos and Baise (2010)
#   PGA1100 = median PGA when Vs30 = 1100 m/s
#   arb = "1" if the standard deviation should be calculated for the arbitrary
#         horizontal component; 0 if the standard deviation should be calculated
#         for the geometric mean of LnY
#   epsilon = number of standard deviations to be considered in the calculations
#   T = Spectral period, sec (0 for PGA; -1 for PGV; -2 for PGD)

# OUTPUT PARAMETERS (from Sa function):
#   Sa = Spectral acceleration (g)
#   PGA = Peak ground acceleration (g); calculated by evaluating Sa at T = 0
#   PGV = Peak ground velocity (cm/sec); calculated by evaluating Sa at T = -1
#   PGD = Peak ground displacement (cm); calculated by evaluating Sa at T = -2




# 1. MODEL COEFFICIENTS

# 1a. Periods with defined coefficients (PGA is 0; PGV is -1; PGD is -2)
periods.cb <- function(positive = FALSE) {

  # Return list of periods excluding PGA, PGV, and PGD
  if(positive){
    T <- c(0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40,
           0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0)

  # Return list of periods including PGA and PGV
  }else{
    T <- c(0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40,
           0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0, 0.0, -1.0, -2.0)
  }
  return(T)
}


# 1b. Period-independent coefficients
#     Bottom of Table 2 in Campbell and Bozorgnia (2008)
coefs.cb <- function(){
  c("c = 1.88", "n = 1.18")
}


# 1c. Period-dependent coefficients for median ground motion term
#     Table 2 in Campbell and Bozorgnia (2008)

c0.cb <- function(T) {
  Period.list <- periods.cb()
  c0.list <- c(-1.715, -1.68, -1.552, -1.209, -0.657, -0.314, -0.133, -0.486,
               -0.89, -1.171, -1.466, -2.569, -4.844, -6.406, -8.692, -9.701,
               -10.556, -11.212, -11.684, -12.505, -13.087, -1.715, 0.954, -5.270)
  c0.list[match(T, Period.list)]
}

c1.cb <- function(T) {
  Period.list <- periods.cb()
  c1.list <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.656, 0.972,
               1.196, 1.513, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 0.5, 0.696, 1.600)
  c1.list[match(T, Period.list)]
}

c2.cb <- function(T) {
  Period.list <- periods.cb()
  c2.list <- c(-0.53, -0.53, -0.53, -0.53, -0.53, -0.53, -0.53, -0.446, -0.362,
               -0.294, -0.186, -0.304, -0.578, -0.772, -1.046, -0.978, -0.638,
               -0.316, -0.07, -0.07, -0.07, -0.53, -0.309, -0.070)
  c2.list[match(T, Period.list)]
}

c3.cb <- function(T) {
  Period.list <- periods.cb()
  c3.list <- c(-0.262, -0.262, -0.262, -0.267, -0.302, -0.324, -0.339, -0.398, -0.458,
               -0.511, -0.592, -0.536, -0.406, -0.314, -0.185, -0.236, -0.491, -0.77,
               -0.986, -0.656, -0.422, -0.262, -0.019, 0.000)
  c3.list[match(T, Period.list)]
}

c4.cb <- function(T) {
  Period.list <- periods.cb()
  c4.list <- c(-2.118, -2.123, -2.145, -2.199, -2.277, -2.318, -2.309, -2.22, -2.146,
               -2.095, -2.066, -2.041, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2.118,
               -2.016, -2.000)
  c4.list[match(T, Period.list)]
}

c5.cb <- function(T) {
  Period.list <- periods.cb()
  c5.list <- c(0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17,
               0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17)
  c5.list[match(T, Period.list)]
}

c6.cb <- function(T) {
  Period.list <- periods.cb()
  c6.list <- c(5.6, 5.6, 5.6, 5.74, 7.09, 8.05, 8.79, 7.6, 6.58, 6.04, 5.3, 4.73, 4,
               4, 4, 4, 4, 4, 4, 4, 4, 5.6, 4, 4)
  c6.list[match(T, Period.list)]
}

c7.cb <- function(T) {
  Period.list <- periods.cb()
  c7.list <- c(0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28,
               0.28, 0.255, 0.161, 0.094, 0, 0, 0, 0, 0, 0.28, 0.245, 0)
  c7.list[match(T, Period.list)]
}

c8.cb <- function(T) {
  Period.list <- periods.cb()
  c8.list <- c(-0.12, -0.12, -0.12, -0.12, -0.12, -0.099, -0.048, -0.012, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, -0.12, 0, 0)
  c8.list[match(T, Period.list)]
}

c9.cb <- function(T) {
  Period.list <- periods.cb()
  c9.list <- c(0.49, 0.49, 0.49, 0.49, 0.49, 0.49, 0.49, 0.49, 0.49, 0.49, 0.49, 0.49,
               0.49, 0.49, 0.49, 0.371, 0.154, 0, 0, 0, 0, 0.49, 0.358, 0)
  c9.list[match(T, Period.list)]
}

c10.cb <- function(T) {
  Period.list <- periods.cb()
  c10.list <- c(1.058, 1.102, 1.174, 1.272, 1.438, 1.604, 1.928, 2.194, 2.351, 2.46,
                2.587, 2.544, 2.133, 1.571, 0.406, -0.456, -0.82, -0.82, -0.82, -0.82,
                -0.82, 1.058, 1.694, -0.82)
  c10.list[match(T, Period.list)]
}

c11.cb <- function(T) {
  Period.list <- periods.cb()
  c11.list <- c(0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04,
                0.077, 0.15, 0.253, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.04, 0.092, 0.3)
  c11.list[match(T, Period.list)]
}

c12.cb <- function(T) {
  Period.list <- periods.cb()
  c12.list <- c(0.61, 0.61, 0.61, 0.61, 0.61, 0.61, 0.61, 0.61, 0.7, 0.75, 0.85, 0.883,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 0.61, 1, 1)
  c12.list[match(T, Period.list)]
}

k1.cb <- function(T) {
  Period.list <- periods.cb()
  k1.list <- c(865, 865, 908, 1054, 1086, 1032, 878, 748, 654, 587, 503, 457, 410, 400,
               400, 400, 400, 400, 400, 400, 400, 865, 400, 400)
  k1.list[match(T, Period.list)]
}

k2.cb <- function(T) {
  Period.list <- periods.cb()
  k2.list <- c(-1.186, -1.219, -1.273, -1.346, -1.471, -1.624, -1.931, -2.188, -2.381,
               -2.518, -2.657, -2.669, -2.401, -1.955, -1.025, -0.299, 0, 0, 0, 0, 0,
               -1.186, -1.955, 0)
  k2.list[match(T, Period.list)]
}

k3.cb <- function(T) {
  Period.list <- periods.cb()
  k3.list <- c(1.839, 1.84, 1.841, 1.843, 1.845, 1.847, 1.852, 1.856, 1.861, 1.865,
               1.874, 1.883, 1.906, 1.929, 1.974, 2.019, 2.11, 2.2, 2.291, 2.517,
               2.744, 1.839, 1.929, 2.744)
  k3.list[match(T, Period.list)]
}


# 1d. Coefficients for standard deviation term
#     Table 3 in Campbell and Bozorgnia (2008)

Sigma.LnY.cb <- function(T) {
  Period.list <- periods.cb()
  Sigma.LnY.list <- c(0.478, 0.480, 0.489, 0.510, 0.520, 0.531, 0.532, 0.534, 0.534,
                      0.544, 0.541, 0.550, 0.568, 0.568, 0.564, 0.571, 0.558, 0.576,
                      0.601, 0.628, 0.667, 0.478, 0.484, 0.667)
  Sigma.LnY.list[match(T, Period.list)]
}

Tau.LnY.cb <- function(T) {
  Period.list <- periods.cb()
  Tau.LnY.list <- c(0.219, 0.219, 0.235, 0.258, 0.292, 0.286, 0.280, 0.249, 0.240,
                    0.215, 0.217, 0.214, 0.227, 0.255, 0.296, 0.296, 0.326, 0.297,
                    0.359, 0.428, 0.485, 0.219, 0.203, 0.485)
  Tau.LnY.list[match(T, Period.list)]
}

Sigma.C.cb <- function(T) {
  Period.list <- periods.cb()
  Sigma.C.list <- c(0.166, 0.166, 0.165, 0.162, 0.158, 0.170, 0.180, 0.186, 0.191,
                    0.198, 0.206, 0.208, 0.221, 0.225, 0.222, 0.226, 0.229, 0.237,
                    0.237, 0.271, 0.290, 0.166, 0.190, 0.290)
  Sigma.C.list[match(T, Period.list)]
}

rho.cb <- function(T) {
  Period.list <- periods.cb()
  rho.list <- c(1.000, 0.999, 0.989, 0.963, 0.922, 0.898, 0.890, 0.871, 0.852,
                0.831, 0.785, 0.735, 0.628, 0.534, 0.411, 0.331, 0.289, 0.261,
                0.200, 0.174, 0.174, 1.000, 0.691, 0.174)
  rho.list[match(T, Period.list)]
}




# 2. NECESSARY FUNCTIONS FOR CALCULATING MEDIAN GROUND MOTION

# 2a. Magnitude Term (Eqn 2)
Fmag.cb <- function(M, T) {
  if(M <= 5.5) {
    Fmag <- c0.cb(T) + c1.cb(T)*M
  } else {
    if(M > 5.5 & M <= 6.5){
      Fmag <- c0.cb(T) + c1.cb(T)*M + c2.cb(T)*(M - 5.5)
    } else {
      if(M > 6.5){
        Fmag <- c0.cb(T) + c1.cb(T)*M + c2.cb(T)*(M - 5.5) + c3.cb(T)*(M - 6.5)
      }
    }
  }
  Fmag
}


# 2b. Distance Term (Eqn 3)
Fdist.cb <- function(M, Rrup, T) {
  return((c4.cb(T) + c5.cb(T)*M)*log(sqrt(Rrup^2 + c6.cb(T)^2)))
}


# 2c. Fault Mechanism Term
Fflt.cb <- function(Ztor, Frv, Fnm, T) {

  # Calculate Fflt.z (Eqn 5)
  if(Ztor < 1) {
    Fflt.z <- Ztor
  } else {
    if(Ztor >= 1) {
      Fflt.z <- 1
    }
  }

  # Return Fflt (Eqn 4)
  return(c7.cb(T)*Frv*Fflt.z + c8.cb(T)*Fnm)
}


# 2d. Hanging-Wall Term
Fhng.cb <- function(M, Rrup, Rjb, Ztor, dip, T) {

  # Calculate Fhng.R (Rupture Distance), Eqn 7
  if(Rjb == 0) {
    Fhng.R <- 1
  } else {
    if(Rjb > 0 & Ztor < 1) {
      Fhng.R <- (max(Rrup, sqrt(Rjb^2 + 1)) - Rjb) / max(Rrup, sqrt(Rjb^2 + 1))
    } else {
      if(Rjb > 0 & Ztor >= 1) {
        Fhng.R <- (Rrup - Rjb)/Rrup
      }
    }
  }

  # Calculate Fhng.M (Magnitude), Eqn 8
  if(M <= 6) {
    Fhng.M <- 0
  } else {
    if(M > 6 & M < 6.5) {
      Fhng.M <- 2*(M - 6)
    } else {
      if(M >= 6.5) {
        Fhng.M <- 1
      }
    }
  }

  # Calculate Fhng.Z (Depth to Top of Rupture), Eqn 9
  if(Ztor >= 20) {
    Fhng.Z <- 0
  } else {
    if(Ztor >= 0 & Ztor < 20) {
      Fhng.Z <- (20 - Ztor)/20
    }
  }

  # Calculate Fhng.dip (Dip Angle), Eqn 10
  if(dip <= 70) {
    Fhng.dip <- 1
  } else {
    if(dip > 70) {
      Fhng.dip <- (90 - dip)/20
    }
  }

  # Return Fhng, Eqn 6
  return( c9.cb(T) * Fhng.R * Fhng.M * Fhng.Z * Fhng.dip )
}


# 2e. Shallow Site Response Term (Eqn 11)
Fsite.cb <- function(Vs30, A1100, T) {

  # Load period-independent coefficients
  coefs.list <- coefs.cb()
  for(i in 1:length(coefs.list))
    eval(parse(text = coefs.list[i]))

  # Calculate site response term
  if(Vs30 < k1.cb(T)){
    Fsite <- c10.cb(T)*log(Vs30/k1.cb(T)) +
      k2.cb(T)*(log(A1100 + c*(Vs30/k1.cb(T))^n) - log(A1100 + c))
  } else {
    if(Vs30 >= k1.cb(T) & Vs30 < 1100){
      Fsite <- (c10.cb(T) + k2.cb(T)*n)*log(Vs30/k1.cb(T))
    } else {
      if(Vs30 >= 1100){
        Fsite <- (c10.cb(T) + k2.cb(T)*n)*log(1100/k1.cb(T))
      }
    }
  }
  Fsite
}


# 2f. Basin Response Term (Eqn 12)
Fsed.cb <- function(Z2.5, T) {
  if(Z2.5 < 1) {
    Fsed <- c11.cb(T)*(Z2.5 - 1)
  } else {
    if(Z2.5 >= 1 & Z2.5 <= 3) {
      Fsed <- 0
    } else {
      if(Z2.5 > 3) {
        Fsed <- c12.cb(T)*k3.cb(T)*exp(-0.75)*(1 - exp(-0.25*(Z2.5-3)))
      }
    }
  }
  Fsed
}


# 2g. Calculation of A1100 (median PGA when Vs30 = 1100 m/s)
A1100.cb <- function(M, Rrup, Rjb, Ztor, Frv, Fnm, dip, Z2.5) {
  LnA1100 <- Fmag.cb(M, T=0) + Fdist.cb(M, Rrup, T=0) +
    Fflt.cb(Ztor, Frv, Fnm, T=0) + Fhng.cb(M, Rrup, Rjb, Ztor, dip, T=0) +
      Fsite.cb(Vs30=1100, A1100=0, T=0) + Fsed.cb(Z2.5, T=0)
  return(exp(LnA1100))
}




# 3. NECESSARY FUNCTIONS FOR CALCULATING STANDARD DEVIATION TERM

# 3a. Partial derivative of site response function (Fsite) with respect to A1100;
#     Alpha (Eqn 17)
Alpha.cb <- function(A1100, Vs30, T){

  # Load period-independent coefficients
  coefs.list <- coefs.cb()
  for(i in 1:length(coefs.list))
    eval(parse(text = coefs.list[i]))

  # Calculate alpha
  if(Vs30 < k1.cb(T))
    k2.cb(T)*A1100 * (1/(A1100 + c*(Vs30/k1.cb(T))^n) - 1/(A1100 + c))
  else
    0
}


# 3b. Sigma, intra-event standard deviation
Sigma.cb <- function(A1100, Vs30, T){

  # Define SigmaLnAF, the estimated SD of the logarithm of the site
  # amplification factor, assuming linear site response
  SigmaLnAF <- 0.3

  # Calculate the standard deviation of ground motion at the base of
  # the site profile
    # At the period of interest (T)
    SigmaLnYb <- sqrt(Sigma.LnY.cb(T)^2 - SigmaLnAF^2)
    # For PGA
    SigmaLnAb <- sqrt(Sigma.LnY.cb(0)^2 - SigmaLnAF^2)

  # Calculate alpha, the partial derivative of the site response function
  # with respect to A1100 (Eqn 17)
  alpha <- Alpha.cb(A1100, Vs30, T)

  # Calculate Sigma, the intra-event standard deviation (Eqn 15)
  sqrt(SigmaLnYb^2 + SigmaLnAF^2 + (alpha^2)*(SigmaLnAb^2) +
       2 * alpha * rho.cb(T) * SigmaLnYb * SigmaLnAb)
}


# 3c. Tau, inter-event standard deviation (Eqn 14)
Tau.cb <- function(T) {
  Tau.LnY.cb(T)
}


# 3d. SigmaTot, total standard deviation of the geometric mean
#     horizontal component of LnY
SigmaTot.cb <- function(A1100, Vs30, T){

  # Intra-event standard deviation (Eqn 15)
  Sigma <- Sigma.cb(A1100, Vs30, T)

  # Inter-event standard deviation (Eqn 14)
  Tau <- Tau.cb(T)

  # Total standard deviation (Eqn 16)
  sqrt(Sigma^2 + Tau^2)
}


# 3e. SigmaArb, total standard deviation of the arbitrary
#     horizontal component of LnY
SigmaArb.cb <- function(A1100, Vs30, T){

  # Total standard deviation of geometric mean
  SigmaTot <- SigmaTot.cb(A1100, Vs30, T)

  # Contribution of component-to-component variability
  SigmaC <- Sigma.C.cb(T)

  # Standard deviation of arbitrary horizontal component (Eqn 18)
  sqrt(SigmaTot^2 + SigmaC^2)
  
}




# 4. MEDIAN GROUND MOTION CALCULATION OF Sa, PGA, and PGV
SaMedian.cb <- function(M, Rjb, Rrup, Ztor, Frv, Fnm, dip, Vs30, Z2.5, T){

  # Calculate A1100, the median PGA when Vs30 = 1100
  A1100 <- A1100.cb(M, Rrup, Rjb, Ztor, Frv, Fnm, dip, Z2.5)

  # Calculate spectral acceleration (Eqn 1)
  Sa <- (exp(Fmag.cb(M, T) + Fdist.cb(M, Rrup, T) + Fflt.cb(Ztor, Frv, Fnm, T) +
             Fhng.cb(M, Rrup, Rjb, Ztor, dip, T) + Fsite.cb(Vs30, A1100, T) +
             Fsed.cb(Z2.5, T)))
    
  # Check for PSA < PGA at short periods (note on page 147
  # of Campbell and Bozorgnia (2008))
  if(T <= 0.25){
    PGA <- (exp(Fmag.cb(M, T=0) + Fdist.cb(M, Rrup, T=0) +
                Fflt.cb(Ztor, Frv, Fnm, T=0) +
                Fhng.cb(M, Rrup, Rjb, Ztor, dip, T=0) +
                Fsite.cb(Vs30, A1100, T=0) + Fsed.cb(Z2.5, T=0)))
    if(Sa < PGA)
      Sa <- PGA
  }
  
  return(Sa)
}




# 5. FINAL FUNCTION FOR CB08 GROUND MOTION CALCULATIONS
Sa.cb <- function(M, Rjb, Rrup = NA, rake = NA, Frv = NA, Fnm = NA, dip = NA,
                  W = NA, Ztor = NA, Vs30, Z1.0 = NA, Z1.5 = NA, Z2.5 = NA, Zhyp = NA,
                  Fhw = NA, azimuth = NA, arb = 0, epsilon, T){

  # If T is a vector, perform calculation for each of the elements
  if(length(T) > 1) {
    return(sapply(T, Sa.cb, M = M, Rjb = Rjb, Rrup = Rrup, rake = rake,
                  Frv = Frv, Fnm = Fnm, dip = dip, W = W, Ztor = Ztor,
                  Vs30 = Vs30, Z1.0 = Z1.0, Z1.5 = Z1.5, Z2.5 = Z2.5, Zhyp = Zhyp,
                  Fhw = Fhw, azimuth = azimuth, arb = arb, epsilon = epsilon))

  # Perform calculation for single value of T:
  } else {


    # 5A.  CHECK INPUT PARAMETERS

    # Check mandatory input parameters
    if(is.na(M) == TRUE | M < 0)
      stop("M must be a positive number")
    if(is.na(Rjb) == TRUE | Rjb < 0)
      stop("Rjb must be a non-negative number")
    if(is.na(Vs30) == TRUE | Vs30 < 0)
      stop("Vs30 must be a positive number")
    if(is.na(epsilon) == TRUE)
      stop("epsilon must be numeric")
    if(is.na(arb) == TRUE | !(arb == 1 | arb == 0))
      stop("arb must be either 1 (to output the standard deviation of the arbitrary horizontal component of ground motion)",
           "\n", "or 0 (to output the standard deviation of the geometric mean horizontal component of ground motion, the default)")
    
    # Check style of faulting parameters
    if(is.na(rake) == TRUE & (is.na(Frv) == TRUE | is.na(Fnm) == TRUE))
      stop("either (1) the rake angle, or (2) both Frv and Fnm must be specified")
    if(is.na(rake) == FALSE & (is.na(Frv) == FALSE | is.na(Fnm) == FALSE) &
       (is.na(Frv) == TRUE | is.na(Fnm) == TRUE))
      stop("either (1) the rake angle, or (2) both Frv and Fnm must be specified")  
    if(is.na(rake) == FALSE & abs(rake) > 180)
      stop("rake angle must be between -180 and 180, inclusive")
    if(is.na(Frv) == FALSE & is.na(Fnm) == FALSE){
      if(Frv == 1 & Fnm == 1)
        stop("either Frv or Fnm may be equal to 1, but both flags cannot be equal to 1")
      if(!(Frv == 1 & Fnm == 0) & !(Frv == 0 & Fnm == 1) & !(Frv == 0 & Fnm == 0))
        stop("Frv must be 1 for reverse faulting (30 <= rake <= 150) and 0 otherwise.",
             "\n", "Fnm must be 1 for normal faulting (-150 <= rake <= -30) and 0 otherwise.")
    }
    # Ensure consistency between rake, Frv, and Fnm
    if(is.na(rake) == FALSE & is.na(Frv) == FALSE & is.na(Fnm) == FALSE){
      if(rake >= 30 & rake <= 150 & !(Frv == 1 & Fnm == 0))
        stop("Inconsistency between rake and style-of-faulting flag variables:", "\n",
             "Frv = 1 and Fnm = 0 for reverse faulting (30 <= rake <= 150)")
      if(rake >= -150 & rake <= -30 & !(Frv == 0 & Fnm == 1))
        stop("Inconsistency between rake and style-of-faulting flag variables:", "\n",
             "Frv = 0 and Fnm = 1 for normal faulting (-120 <= rake <= -60)")
      if((abs(rake) < 30 | abs(rake) > 150) & !(Frv == 0 & Fnm == 0))
        stop("Inconsistency between rake and style-of-faulting flag variables:", "\n",
             "Frv = 0 and Fnm = 0 for strike-slip faulting.")
    }

    # Check hanging wall parameters (only necessary when Rrup is not specified, since
    # Rrup is calculated from Rx, which is calculated from Rjb)
    if(is.na(Rrup) == TRUE | Rrup < 0){
      if(is.na(azimuth) == FALSE & abs(azimuth) > 180)
        stop("the source-to-site azimuth must be between -180 and 180, inclusive")
      if(is.na(Fhw) == FALSE & !(Fhw == 1 | Fhw == 0))
        stop("Fhw must be either 1 (for sites on the hanging wall side of the fault)", "\n",
             "or 0 (for sites on the footwall side of the fault)")
      # Ensure consistency between azimuth and Fhw
      if(is.na(azimuth) == FALSE & is.na(Fhw) == FALSE){
        if(azimuth < 0 & Fhw == 1)
          stop("Inconsistency between azimuth and Fhw. Fhw must be 0 when azimuth < 0.")
        if(azimuth > 0 & Fhw == 0)
          stop("Inconsistency between azimuth and Fhw. Fhw must be 1 when azimuth > 0.")
      }
    }


    # 5B. OBTAIN ESTIMATES OF UNSPECIFIED INPUT PARAMETERS
    
    # Assign generic rake angle if Fhw and Fnm are specified
    if(is.na(rake) == TRUE){
      if(Frv == 1 & Fnm == 0){  # Reverse faulting
        rake <- 90
      } else{
        if(Frv == 0 & Fnm == 1){  # Normal faulting
          rake <- -90
        } else{
          if(Frv == 0 & Fnm == 0){  # Strike-slip faulting
            rake <- 180
          }
        }
      }
    }

    # Convert rake to fault type if rake is provided
    if(is.na(Frv) == TRUE & is.na(Fnm) == TRUE){
      # Frv
      if(rake >= 30 & rake <= 150)
        Frv <- 1
      else
        Frv <- 0
      # Fnm
      if(rake >= -150 & rake <= -30)
        Fnm <- 1
      else
        Fnm <- 0
    }

    # Dip angle
    if(is.na(dip) == TRUE | dip < 0)
      dip <- dip.calc(rake)

    # Depth to top of rupture, Ztor
    if(is.na(Ztor) == TRUE | Ztor < 0){
      # Down-dip rupture width, W (used for calculating Ztor)
      if(is.na(W) == TRUE | W < 0)
        W <- W.calc(M, rake)
      Ztor <- Ztor.calc(Zhyp, W, dip, M, rake)
    }

    # Calculate Rrup from Rjb (using distance equations from Kaklamanos and
    # Baise (2010) if Rrup is not provided)
    if(is.na(Rrup) == TRUE | Rrup < 0){
      # Down-dip rupture width, W (used for calculating Rrup)
      if(is.na(W) == TRUE | W < 0)
        W <- W.calc(M, rake)     
      # Azimuth angle (used for calculating Rrup)
      if(is.na(azimuth) == TRUE){
        if(is.na(Fhw) == FALSE){
          # If hanging wall site, assume azimuth = 50 deg
          if(Fhw == 1)
            azimuth <- 50
          # If footwall site, assume azimuth = -50 deg  
          else
            azimuth <- -50         
        # If Fhw is not provided, assume azimuth = -50 deg (footwall site)
        # Rrup is calculated as sqrt(Ztor^2 + Rjb^2)
        } else{
          azimuth <- -50
        }         
      }
      # Site coordinate, Rx (used for calculating Rrup)
      Rx <- Rx.calc(Rjb, Ztor, W, dip, azimuth, Rrup)
      # Rupture distance, Rrup
      Rrup <- Rrup.calc(Rx, Ztor, W, dip, azimuth, Rjb)
    }


    # Depth parameter, Z2.5
    if(is.na(Z2.5) == TRUE){
      # Calculate from Z1.5 if provided
      # (Eqn 6.4 in Campbell and Bozorgnia (2007); final report to PEER)
      if(is.na(Z1.5) == FALSE & Z1.5 > 0){
        Z2.5 <- 0.636 + 0.001549*Z1.5
        # Calculate from Z1.0 if provided
        # (Eqn 6.3 in Campbell and Bozorgnia (2007); final report to PEER)
      } else{
        if(is.na(Z1.0) == FALSE & Z1.0 > 0){
          Z2.5 <- 0.519 + 0.003595*Z1.0
          # If neither Z1.0 nor Z1.5 is provided, estimate from
          # Vs30 using the AS08 relation for Z1.0 = f(Vs30)
        } else{
          Z1.0 <- Z1.calc.as(Vs30)
          Z2.5 <- 0.519 + 0.003595*Z1.0
        }
      }
    }

    
    # 5C. CALCULATE GROUND MOTION PARAMETER
      
    # Is interpolation necessary?
    interp <- getPeriod(T, "CB08")$interp

    # If interpolation is not necessary, compute Sa
    if(interp == FALSE){
      LnSaMedian <- log(SaMedian.cb(M, Rjb, Rrup, Ztor, Frv, Fnm, dip, Vs30, Z2.5, T))
      A1100 <- A1100.cb(M, Rrup, Rjb, Ztor, Frv, Fnm, dip, Z2.5)
      if(is.na(arb) == FALSE & arb == 1)   # Uncertainty of arbitrary horiz. component
        epsilon.SD <- epsilon * SigmaArb.cb(A1100, Vs30, T)
      else    # Uncertainty of geometric mean component
        epsilon.SD <- epsilon * SigmaTot.cb(A1100, Vs30, T)
      LnSa <- LnSaMedian + epsilon.SD
      return(exp(LnSa))
    } else{

    # If interpolation is necessary, compute Sa
      if(interp == TRUE){
        T1 <- getPeriod(T, "CB08")$lower
        T2 <- getPeriod(T, "CB08")$upper

        # PGA1100
        A1100 <- A1100.cb(M, Rrup, Rjb, Ztor, Frv, Fnm, dip, Z2.5)    
      
        # Calculation for T1
        LnSaMedian.T1 <- log(SaMedian.cb(M, Rjb, Rrup, Ztor, Frv, Fnm, dip, Vs30, Z2.5, T1))
        A1100 <- A1100.cb(M, Rrup, Rjb, Ztor, Frv, Fnm, dip, Z2.5)
        if(is.na(arb) == FALSE & arb == 1)   # Uncertainty of arbitrary horiz. component
          epsilon.SD.T1 <- epsilon * SigmaArb.cb(A1100, Vs30, T1)
        else    # Uncertainty of geometric mean component
          epsilon.SD.T1 <- epsilon * SigmaTot.cb(A1100, Vs30, T1)
        LnSaT1 <- LnSaMedian.T1 + epsilon.SD.T1
        
        # Calculation for T2
        LnSaMedian.T2 <- log(SaMedian.cb(M, Rjb, Rrup, Ztor, Frv, Fnm, dip, Vs30, Z2.5, T2))
        A1100 <- A1100.cb(M, Rrup, Rjb, Ztor, Frv, Fnm, dip, Z2.5)
        if(is.na(arb) == FALSE & arb == 1)   # Uncertainty of arbitrary horiz. component
          epsilon.SD.T2 <- epsilon * SigmaArb.cb(A1100, Vs30, T2)
        else    # Uncertainty of geometric mean component
          epsilon.SD.T2 <- epsilon * SigmaTot.cb(A1100, Vs30, T2)
        LnSaT2 <- LnSaMedian.T2 + epsilon.SD.T2
      
        # Interpolated value
        LnSa <- interpolate(log(T), log(T1), log(T2), LnSaT1, LnSaT2)
        return(exp(LnSa))
      }
    }
  }
}

