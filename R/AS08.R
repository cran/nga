# ABRAHAMSON & SILVA NGA MODEL
# Abrahamson, N., and W. Silva (2008). Summary of the Abrahamson & Silva NGA Ground-Motion Relations.
# Earthquake Spectra, Vol. 24, pp. 67-97.


# OUTLINE OF CODE
# 1. Model Coefficients
#    a. Periods with defined coefficients
#    b. Period-independent coefficients
#    c. Period-dependent coefficients for median ground motion term
#    d. Coefficients for standard deviation term
# 2. Necessary Functions for Calculating Median Ground Motion
#    a. Base Model
#    b. Site Response Model
#    c. Hanging Wall Model
#    d. Depth to Top of Rupture Model
#    e. Large Distance Model
#    f. Soil Depth Model
#    g. Calculation of PGA1100
# 3. Necessary Functions for Calculating Standard Deviation Term
#    a. Partial derivative of site response function (f5) with respect to PGA1100
#    b. Sigma0
#    c. Tau0
#    d. Sigma, intra-event standard deviation
#    e. Tau, inter-event standard deviation
#    f. SigmaTot, total standard deviation
# 4. Median Ground Motion Calculation of Sa, PGA, and PGV
# 5. Final Function for AS08 Ground Motion Calculations
#    a. Check input parameters
#    b. Obtain estimates of unspecified input parameters
#    c. Calculate ground motion parameter
#    **** This is the function that users will primarily use for ground motion calculations ****


# INPUT PARAMETERS FOR THE FUNCTIONS BELOW:
#   M = Moment magnitude
#   Rjb = Joyner-Boore distance (km)
#   Vs30 = Time-averaged shear wave velocity over 30 m subsurface depth  (m/sec)
#   VsFlag = Flag variable indicating how Vs30 is obtained
#            (1 if measured, 0 if estimated/inferred)
#   epsilon = number of standard deviations to be considered in the calculations
#   T = Spectral period, sec (0 for PGA; -1 for PGV)
#   Rrup = Rupture distance (km)
#   Rx = Site coordinate (km)
#   dip = Fault dip (deg)
#   W = Down-dip rupture width (km)
#   Ztor = Depth to top of rupture (km)
#   Z1.0 = Depth to Vs = 1.0 km/sec  (m)
#   rake = Rake angle of fault movement (deg)
#   Frv = Reverse style-of-faulting flag (1 for reverse faulting,
#         0 otherwise); calculated from rake
#   Fnm = Normal style-of-faulting flag (1 for normal faulting,
#         0 otherwise); calculated from rake
#   Fhw = Hanging wall flag (1 for site on hanging wall side of fault,
#         0 otherwise); calculated from azimuth or Rx (if provided)
#   azimuth = source-to-site azimuth (deg); see Kaklamanos et al. (2011)
#   Zhyp = hypocentral depth (km)
#   Fas = Aftershock flag (1 for aftershocks, 0 for mainshocks)
#   PGA1100 = median PGA when Vs30 = 1100 m/s


# OUTPUT PARAMETERS (from Sa function):
#   Sa = Spectral acceleration (g)
#   PGA = Peak ground acceleration (g); calculated by evaluating Sa at T = 0
#   PGV = Peak ground velocity (cm/sec); calculated by evaluating Sa at T = -1




# 1. MODEL COEFFICIENTS

# 1a. Periods with defined coefficients (PGA is 0; PGV is -1)
periods.as <- function(positive = FALSE) {
  
  # Return list of periods excluding PGA and PGV
  if(positive){
    T <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.10, 0.15, 0.20, 0.25,
           0.30, 0.40, 0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0)

  # Return list of periods including PGA and PGV  
  }else{
    T <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.10, 0.15, 0.20, 0.25,
           0.30, 0.40, 0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0, 0, -1)
  }
  return(T)
}


# 1b. Period-independent coefficients
#     Table 4 in Abrahamson and Silva (2008)
coefs.as <- function() {
  c("c = 1.88", "c1 = 6.75", "c2 = 50", "c4 = 4.5", "n = 1.18")
}


# 1c. Period-dependent coefficients for median ground motion terms
#     Tables 5a and 5b in Abrahamson and Silva (2008)

Vlin.as <- function(T) {
  Period.list <- periods.as()
  Vlin.list <- c(865.1, 865.1, 907.8, 994.5, 1053.5, 1085.7, 1032.5, 877.6, 748.2,
                 654.3, 587.1, 503.0, 456.6, 410.5, 400.0, 400.0, 400.0, 400.0, 400.0,
                 400.0, 400.0, 400.0, 865.1, 400.0)  
  Vlin.list[match(T, Period.list)]
}

b.as <- function(T) {
  Period.list <- periods.as()
  b.list <- c(-1.186, -1.219, -1.273, -1.308, -1.346, -1.471, -1.624, -1.931, -2.188,
              -2.381, -2.518, -2.657, -2.669, -2.401, -1.955, -1.025, -0.299, 0.000,
              0.000, 0.000, 0.000, 0.000, -1.186, -1.955)  
  b.list[match(T, Period.list)]
}

a1.as <- function(T) {
  Period.list <- periods.as()
  a1.list <- c(0.811, 0.855, 0.962, 1.037, 1.133, 1.375, 1.563, 1.716, 1.687, 1.646,
               1.601, 1.511, 1.397, 1.137, 0.915, 0.510, 0.192, -0.280, -0.639,
               -0.936, -1.527, -1.993, 0.804, 5.7578)  
  a1.list[match(T, Period.list)]
}

a2.as <- function(T) {
  Period.list <- periods.as()
  a2.list <- c(-0.9679, -0.9774, -1.0024, -1.0289, -1.0508, -1.0810, -1.0833, -1.0357,
               -0.9700, -0.9202, -0.8974, -0.8677, -0.8475, -0.8206, -0.8088, -0.7995,
               -0.7960, -0.7960, -0.7960, -0.7960, -0.7960, -0.7960, -0.9679, -0.9046)  
  a2.list[match(T, Period.list)]
}

a3.as <- function(T) {
  Period.list <- periods.as()
  a3.list <- c(0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265,
               0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265,
               0.265, 0.265, 0.265, 0.265)  
  a3.list[match(T, Period.list)]
}

a4.as <- function(T) {
  Period.list <- periods.as()
  a4.list <- c(-0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231,
             -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231,
             -0.231, -0.231, -0.231, -0.231, -0.231, -0.231)  
  a4.list[match(T, Period.list)]
}

a5.as <- function(T) {
  Period.list <- periods.as()
  a5.list<- c(-0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398,
              -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398,
              -0.398, -0.398, -0.398, -0.398, -0.398, -0.398)  
  a5.list[match(T, Period.list)]
}

a8.as <- function(T) {
  Period.list <- periods.as()
  a8.list <- c(-0.0372, -0.0372, -0.0372, -0.0315, -0.0271, -0.0191, -0.0166, -0.0254,
               -0.0396, -0.0539, -0.0656, -0.0807, -0.0924, -0.1137, -0.1289, -0.1534,
               -0.1708, -0.1954, -0.2128, -0.2263, -0.2509, -0.2683, -0.0372, -0.1200)  
  a8.list[match(T, Period.list)]
}

a10.as <- function(T) {
  Period.list <- periods.as()
  a10.list <- c(0.9445, 0.9834, 1.0471, 1.0884, 1.1333, 1.2808, 1.4613, 1.8071, 2.0773,
                2.2794, 2.4201, 2.5510, 2.5395, 2.1493, 1.5705, 0.3991, -0.6072,
                -0.9600, -0.9600, -0.9208, -0.7700, -0.6630, 0.9445, 1.5390)  
  a10.list[match(T, Period.list)]
}

a12.as <- function(T) {
  Period.list <- periods.as()
  a12.list <- c(0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0181, 0.0309,
                0.0409, 0.0491, 0.0619, 0.0719, 0.0800, 0.0800, 0.0800, 0.0800, 0.0800,
                0.0800, 0.0800, 0.0800, 0.0800, 0.0000, 0.0800)  
  a12.list[match(T, Period.list)]
}

a13.as <- function(T) {
  Period.list <- periods.as()
  a13.list <- c(-0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600,
                -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600,
                -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600)  
  a13.list[match(T, Period.list)]
}

a14.as <- function(T) {
  Period.list <- periods.as()
  a14.list <- c(1.0800, 1.0800, 1.1331, 1.1708, 1.2000, 1.2000, 1.2000, 1.1683, 1.1274,
                1.0956, 1.0697, 1.0288, 0.9971, 0.9395, 0.8985, 0.8409, 0.8000,
                0.4793, 0.2518, 0.0754, 0.0000, 0.0000, 1.0800, 0.7000)  
  a14.list[match(T, Period.list)]
}

a15.as <- function(T) {
  Period.list <- periods.as()
  a15.list <- c(-0.3500, -0.3500, -0.3500, -0.3500, -0.3500, -0.3500, -0.3500, -0.3500,
                -0.3500, -0.3500, -0.3500, -0.3500, -0.3191, -0.2629, -0.2230, -0.1668,
                -0.1270, -0.0708, -0.0309, 0.0000, 0.0000, 0.0000, -0.3500, -0.3900)  
  a15.list[match(T, Period.list)]
}

a16.as <- function(T) {
  Period.list <- periods.as()
  a16.list <- c(0.9000, 0.9000, 0.9000, 0.9000, 0.9000, 0.9000, 0.9000, 0.9000, 0.9000,
                0.9000, 0.9000, 0.8423, 0.7458, 0.5704, 0.4460, 0.2707, 0.1463,
                -0.0291, -0.1535, -0.2500, -0.2500, -0.2500, 0.9000, 0.6300)  
  a16.list[match(T, Period.list)]
}

a18.as <- function(T) {
  Period.list <- periods.as()
  a18.list <- c(-0.0067, -0.0067, -0.0067, -0.0067, -0.0076, -0.0093, -0.0093,
                -0.0093, -0.0083, -0.0069, -0.0057, -0.0039, -0.0025, 0.0000, 0.0000,
                0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, -0.0067, 0.0000)  
  a18.list[match(T, Period.list)]
}


# 1d. Coefficients for standard deviation term
#     Table 6 in Abrahamson and Silva (2008)

s1.est.as <- function(T) {   # Estimated Vs30
  Period.list <- periods.as()
  s1.est.list <- c(0.590, 0.590, 0.605, 0.615, 0.623, 0.630, 0.630, 0.630, 0.630, 0.630,
                   0.630, 0.630, 0.630, 0.630, 0.630, 0.615, 0.604, 0.589, 0.578, 0.570,
                   0.611, 0.640, 0.590, 0.590)
  s1.est.list[match(T, Period.list)]
}

s2.est.as <- function(T) {   # Estimated Vs30
  Period.list <- periods.as()
  s2.est.list <- c(0.470, 0.470, 0.478, 0.483, 0.488, 0.495, 0.501, 0.509, 0.514, 0.518,
                   0.522, 0.527, 0.532, 0.539, 0.545, 0.552, 0.558, 0.565, 0.570, 0.587,
                   0.618, 0.640, 0.470, 0.470)
  s2.est.list[match(T, Period.list)]
}

s1.meas.as <- function(T) {   # Measured Vs30
  Period.list <- periods.as()
  s1.meas.list <- c(0.576, 0.576, 0.591, 0.602, 0.610, 0.617, 0.617, 0.616, 0.614, 0.612,
                    0.611, 0.608, 0.606, 0.602, 0.594, 0.566, 0.544, 0.527, 0.515, 0.510,
                    0.572, 0.612, 0.576, 0.576)
  s1.meas.list[match(T, Period.list)]
}

s2.meas.as <- function(T) {   # Measured Vs30
  Period.list <- periods.as()
  s2.meas.list <- c(0.453, 0.453, 0.461, 0.466, 0.471, 0.479, 0.485, 0.491, 0.495, 0.497,
                    0.499, 0.501, 0.504, 0.506, 0.503, 0.497, 0.491, 0.500, 0.505, 0.529,
                    0.579, 0.612, 0.453, 0.453)
  s2.meas.list[match(T, Period.list)]
}

s3.as <- function(T) {
  Period.list <- periods.as()
  s3.list <- c(0.420, 0.420, 0.462, 0.492, 0.515, 0.550, 0.550, 0.550, 0.520, 0.497,
               0.479, 0.449, 0.426, 0.385, 0.350, 0.350, 0.350, 0.350, 0.350, 0.350,
               0.350, 0.350, 0.470, 0.420)
  s3.list[match(T, Period.list)]
}

s4.as <- function(T) {
  Period.list <- periods.as()
  s4.list <- c(0.300, 0.300, 0.305, 0.309, 0.312, 0.317, 0.321, 0.326, 0.329, 0.332,
               0.335, 0.338, 0.341, 0.346, 0.350, 0.350, 0.350, 0.350, 0.350, 0.350,
               0.350, 0.350, 0.300, 0.300)
  s4.list[match(T, Period.list)]
}

rho.as <- function(T) {
  Period.list <- periods.as()
  rho.list <- c(1.000, 1.000, 0.991, 0.982, 0.973, 0.952, 0.929, 0.896, 0.874, 0.856,
                0.841, 0.818, 0.783, 0.680, 0.607, 0.504, 0.431, 0.328, 0.255, 0.200,
                0.200, 0.200, 1.000, 0.740)
  rho.list[match(T, Period.list)]
}




# 2. NECESSARY FUNCTIONS FOR CALCULATING MEDIAN GROUND MOTION

# 2a. Base Model
f1.as <- function(M, Rrup, T) {

  # Load period-independent coefficients
  coefs.list <- coefs.as()
  for(i in 1:length(coefs.list))
    eval(parse(text = coefs.list[i]))
  
  # Calculate R (Eqn 3)
  R <- sqrt(Rrup^2 + c4^2)

  # Return f1 (Eqn 2)
  if(M <= c1) {
    f1 <- a1.as(T) + a4.as(T)*(M - c1) + a8.as(T)*(8.5 - M)^2 +
      (a2.as(T) + a3.as(T)*(M - c1))*log(R)
  } else{
    if(M > c1){
      f1 <- a1.as(T) + a5.as(T)*(M - c1) + a8.as(T)*(8.5 - M)^2 +
            (a2.as(T) + a3.as(T)*(M - c1))*log(R)
    }
  }
  f1
}


# 2b. Site Response Model
f5.as <- function(PGA1100, Vs30, T) {

  # Load period-independent coefficients
  coefs.list <- coefs.as()
  for(i in 1:length(coefs.list))
    eval(parse(text = coefs.list[i]))

  # Calculate V1 (Eqn 6)
  if(T == -1) { # PGV
    V1 <- 862
  } else {
    if(T <= 0.5) {
      V1 <- 1500
    } else {
      if(T > 0.5 & T <= 1) {
        V1 <- exp(8 - 0.795*log(T/0.21))
      } else { 
        if(T > 1 & T < 2) {
          V1 <- exp(6.76 - 0.297*log(T))
        } else {
          if(T >= 2)
            V1 <- 700
        }
      }
    }
  }
    
  # Calculate Vs30star (Eqn 5)
  if(Vs30 < V1) {Vs30star <- Vs30} else
  if(Vs30 >= V1) {Vs30star <- V1}
        
  # Return f5 (Eqn 4)
  if(Vs30 < Vlin.as(T)) {
    f5 <- a10.as(T)*log(Vs30star/Vlin.as(T)) - b.as(T)*log(PGA1100 + c) +
            b.as(T)*log(PGA1100 + c*(Vs30star/Vlin.as(T))^n)
  } else {
    if(Vs30 >= Vlin.as(T)) {
      f5 <- (a10.as(T) + b.as(T)*n)*log(Vs30star/Vlin.as(T))
    }
  }
  f5
}


# 2c. Hanging Wall Model
f4.as <- function(Rjb, Rx, dip, Ztor, M, W, T) {

  # Exit if Rx < 0 (site is located on footwall)
  if(Rx < 0)
    return(0)

  else{
    
    # Calculate tapers
    
    # T1 (Eqn 8)
    if(Rjb < 30) {
      T1 <- 1 - Rjb/30
    } else {
      if(Rjb >= 30) {
        T1 <- 0
      }
    }
    
    # T2 (Eqn 9)
    if(Rx > W*cos(dip*pi/180) | dip == 90) {
      T2 <- 1
    } else {
      T2 <- 0.5 + Rx/(2*W*cos(dip*pi/180))
    }
    
    # T3 (Eqn 10)
    if(Rx >= Ztor) {
      T3 <- 1
    } else {
      if(Rx < Ztor) {
        T3 <- Rx/Ztor
      }
    }
    
    # T4 (Eqn 11)
    if(M <= 6) {
      T4 <- 0
    } else { 
      if(M > 6 & M < 7) {
        T4 <- M - 6
      } else {
        if(M >= 7) {
          T4 <- 1
        }
      }
    }
    
    # T5 (Eqn 12 in Abrahamson and Silva (2008), modified by Eqn 5 in
    # Abrahamson and Silva (2009), errata to the AS08 model

    if(dip >= 30) {
      T5 <- 1 - (dip - 30)/60
    } else {
      if(dip < 30) {
        T5 <- 1
      }
    }  
  
    # Return f4 (Eqn 7)
    return(a14.as(T)*T1*T2*T3*T4*T5)
  }
}


# 2d. Depth to Top of Rupture Model (Eqn 13)
f6.as <- function(Ztor, T) {
  if(Ztor < 10) {
    return(a16.as(T)*Ztor/10)
  } else {
    if(Ztor >= 10) {
      return(a16.as(T))
    }
  }
}


# 2e. Large Distance Model
f8.as <- function(Rrup, M, T) {
  
  # Calculate Taper 6 (Eqn 15)
  if(M < 5.5) {
    T6 <- 1
  } else {
    if(M >= 5.5 & M <= 6.5) {
      T6 = 0.5*(6.5 - M) + 0.5
    } else{
      if(M > 6.5) {
        T6 <- 0.5
      }
    }
  }
       
  # Return f8 (Eqn 14)
  if(Rrup < 100) {
    return(0)
  } else {
    if(Rrup >=100) {
      return(a18.as(T)*(Rrup - 100)*T6)
    }
  }
}


# 2f. Soil Depth Model
f10.as <- function(Z1.0, Vs30, T) {

  # Load period-independent coefficients
  coefs.list <- coefs.as()
  for(i in 1:length(coefs.list))
    eval(parse(text = coefs.list[i]))

  # Calculate V1 (Eqn 6)
  if(T == -1) {  # PGV
    V1 <- 862
  } else { 
    if(T <= 0.5) {
      V1 <- 1500
    } else {
      if(T > 0.5 & T <= 1) {
        V1 <- exp(8 - 0.795*log(T/0.21))
      } else {
        if(T > 1 & T < 2) {
          V1 <- exp(6.76 - 0.297*log(T))
        } else {
          if(T >= 2) V1 <- 700
        }
      }
    }
  }
    
  # Calculate Vs30star (Eqn 5)
  if(Vs30 < V1) {
    Vs30star <- Vs30
  } else {
    Vs30star <- V1
  }
        
  # Obtain Ln of Median Z1.0 (Eqn 17)
  if(Vs30 < 180) {
    LnMedianZ1.0 <- 6.745
  } else {
    if(Vs30 >= 180 & Vs30 <= 500) {
      LnMedianZ1.0 <- 6.745 - 1.35*log(Vs30/180)
    } else{
      if(Vs30 > 500) {
        LnMedianZ1.0 <- 5.394 - 4.48*log(Vs30/500)
      }
    }
  }
    
  # Obtain Median Z1.0
  MedianZ1.0 <- exp(LnMedianZ1.0)
    
  # Calculate e2 (Eqn 19)
  if(T == -1) {  # PGV
    e2 <- -0.25 * log(Vs30/1000) * log(1/0.35)
  } else {
    if(T < 0.35 | Vs30 > 1000) {
      e2 <- 0
    } else {
      if(T >= 0.35 & T <= 2) {
        e2 <- -0.25 * log(Vs30/1000) * log(T/0.35)
      } else {
        if(T > 2) {
          e2 <- -0.25 * log(Vs30/1000) * log(2/0.35)
        }
      }
    }
  }
  
  # Calculate a21 (Eqn 18)
  if(Vs30 >= 1000) {
    a21 <- 0
  } else {
    if((a10.as(T) + b.as(T)*n)*log(Vs30star / min(V1, 1000)) +
       e2*log((Z1.0 + c2)/(MedianZ1.0 + c2)) < 0) {
      a21 <- -(a10.as(T) + b.as(T)*n) *
        log(Vs30star / min(V1, 1000))/log((Z1.0 + c2)/(MedianZ1.0 + c2))
    } else
      a21 <- e2
  }
  # Calculate a22 (Eqn 20)
  if(T < 2) {
    a22 <- 0
  } else {
    if(T >= 2) {
      a22 <- 0.0625*(T - 2)
    }
  }
    
  # Return f10 (Eqn 16)
  if(Z1.0 >= 200) {
    f10 <- a21*log((Z1.0 + c2) / (MedianZ1.0 + c2)) + a22*log(Z1.0/200)
  }else{
    if(Z1.0 < 200) {
      f10 <- a21*log((Z1.0 + c2) / (MedianZ1.0 + c2))
    }
  }
  f10
}


# 2g. Calculation of PGA1100 (median PGA when Vs30 = 1100 m/s)
PGA1100.as <- function(M, Rrup, Rjb, Rx, Ztor, Frv, Fnm, Fas, Fhw, dip, W){

  # Define parameters for rock site
  Vs30.rock <- 1100
  Z1.0.rock <- Z1.calc.as(Vs30 = 1100)
  PGA1100.rock <- 0

  # Calculate PGA1100
  LnPGA.rock <- f1.as(M, Rrup, T = 0) + a12.as(T = 0)*Frv +
    a13.as(T = 0) * Fnm + a15.as(T = 0) * Fas +
      f5.as(PGA1100.rock, Vs30.rock, T = 0) +
        Fhw*f4.as(Rjb, Rx, dip, Ztor, M, W, T = 0) +
          f6.as(Ztor, T = 0) + f8.as(Rrup, M, T = 0) +
            f10.as(Z1.0.rock, Vs30.rock, T = 0)
  return(exp(LnPGA.rock))
}




# 3. NECESSARY FUNCTIONS FOR CALCULATING STANDARD DEVIATION TERM

# 3a. Partial derivative of site response function (f5) with respect to PGA1100;
#     Alpha (Eqn 26; corrected in errata per Abrahamson and Silva, 2009)
Alpha.as <- function(PGA1100, Vs30, T){

  # Load period-independent coefficients
  coefs.list <- coefs.as()
  for(i in 1:length(coefs.list))
    eval(parse(text = coefs.list[i]))

  # Calculate alpha  
  if(Vs30 >= Vlin.as(T)) {
    0
  } else{
    if(Vs30 < Vlin.as(T)){
      (-b.as(T)*PGA1100 / (PGA1100 + c)) +
        (b.as(T)*PGA1100 / (PGA1100 + c*(Vs30/Vlin.as(T))^n))
    }
  }
}


# 3b. Sigma0, intra-event standard deviation of the observed
#     ground motion for the linear site response range (Eqn 27)
Sigma0.as <- function(M, VsFlag, T) {
 
  # Obtain coefficients s1 and s2, depending on whether
  # Vs30 is estimated (inferred) or measured
  if(VsFlag == 0){   # Inferred Vs30
    s1 <- s1.est.as(T)
    s2 <- s2.est.as(T)
  } else{
    if(VsFlag == 1){   # Measured Vs30
      s1 <- s1.meas.as(T)
      s2 <- s2.meas.as(T)
    }
  }

  # Calculate Sigma0 (Eqn 27)
  if(M < 5){
    Sigma0 <- s1
  } else{
    if(M >= 5 & M <= 7){
      Sigma0 <- s1 + ((s2-s1)/2)*(M-5)
    } else{
      if(M > 7){
        Sigma0 <- s2
      }
    }
  }
}


# 3c. Tau0, inter-event standard deviation of the observed ground
#     motion for the linear site response range (Eqn 28)
Tau0.as <- function(M, T) {
  if(M < 5){
    Tau0 <- s3.as(T)
  } else{
    if(M >= 5 & M <= 7){
      Tau0 <- s3.as(T) + ((s4.as(T)-s3.as(T))/2)*(M-5)
    } else{
      if(M > 7){
        Tau0 <- s4.as(T)
      }
    }
  }
}


# 3d. Sigma, intra-event standard deviation
Sigma.as <- function(M, PGA1100, Vs30, VsFlag, T){

  # Define SigmaAmp, the intra-event standard deviation of
  # the site amplification factors
  SigmaAmp <- 0.3

  # Calculate SigmaB, the intra-event standard deviation of
  # the input rock motion (Eqn 23)
    # At the period of interest (T)
    SigmaB.T <- sqrt(Sigma0.as(M, VsFlag, T)^2 - SigmaAmp^2)
    # For PGA
    SigmaB.PGA <- sqrt(Sigma0.as(M, VsFlag, 0)^2 - SigmaAmp^2)

  # Calculate alpha, the partial derivative of the site response
  # function (f5) with respect to PGA1100 (Eqn 26)
  alpha <- Alpha.as(PGA1100, Vs30, T)

  # Calculate Sigma, the intra-event standard deviation
  # (Eqn 24; corrected in errata per Abrahamson and Silva, 2009)
  sqrt(SigmaB.T^2 + SigmaAmp^2 + (alpha^2)*(SigmaB.PGA^2) +
       2 * alpha * SigmaB.T * SigmaB.PGA * rho.as(T))
}


# 3e. Tau, inter-event standard deviation
Tau.as <- function(M, PGA1100, Vs30, T){
                 
  # Calculate Tau0, the inter-event standard deviation of the
  # observed ground motion for the linear site response range (Eqn 28)
  Tau0 <- Tau0.as(M, T)

  # Calculate TauB, the inter-event standard deviation of the input rock motion
    # At the period of interest (T)
    TauB.T <- Tau0.as(M, T)
    # For PGA
    TauB.PGA <- Tau0.as(M, 0)

  # Calculate alpha, the partial derivative of the
  # site response function (f5) with respect to PGA1100 (Eqn 26)
  alpha <- Alpha.as(PGA1100, Vs30, T)

  # Calculate Tau, the inter-event standard deviation (Eqn 25)
  sqrt(Tau0^2 + (alpha^2)*(TauB.PGA^2) +
       2 * alpha * TauB.T * TauB.PGA * rho.as(T))
}  


# 3f. SigmaTot, total standard deviation
SigmaTot.as <- function(M, PGA1100, Vs30, VsFlag, T){

  # Intra-event standard deviation
  Sigma <- Sigma.as(M, PGA1100, Vs30, VsFlag, T)

  # Inter-event standard deviation
  Tau <- Tau.as(M, PGA1100, Vs30, T)
  
  # Total standard deviation
  sqrt(Sigma^2 + Tau^2)
}




# 4. MEDIAN GROUND MOTION CALCULATION of Sa, PGA, and PGV
SaMedian.as <- function(M, Rjb, Rrup, Rx, Ztor, Frv, Fnm, Fas, Fhw, dip, Vs30, Z1.0, W, T){
   
  # Calculation of PGA1100, the median PGA when Vs30 = 1100 m/s
  PGA1100 <- PGA1100.as(M, Rrup, Rjb, Rx, Ztor, Frv, Fnm, Fas, Fhw, dip, W)
  Vs30.rock <- 1100
  Z1.0.rock <- Z1.calc.as(Vs30.rock)
    
  # Calculate cutoff period for constant displacement Model (Eqn 21)
  LogTd <- -1.25 + 0.3*M
  Td <- min(10^(LogTd), 10)
    
  # Do not apply constant displacement model for T <= Td
  if(T <= Td){
    LnSa <- f1.as(M, Rrup, T) + a12.as(T) * Frv + a13.as(T) * Fnm +
      a15.as(T)*Fas + f5.as(PGA1100, Vs30, T) +
        Fhw*f4.as(Rjb, Rx, dip, Ztor, M, W, T) + f6.as(Ztor, T) +
          f8.as(Rrup, M, T) + f10.as(Z1.0, Vs30, T)
    Sa <- exp(LnSa)
  } else {
    # Apply constant displacement model for T > Td
    # (modification of Eqn 22 to fix an error in the printed equation)
    
    # Calculate LnSa.rock.Td (Rock Sa at Td for Vs30 = 1100 m/s)
    # Lower bound Td for interpolation
    Td1 <- getPeriod(Td, "AS08")$lower
    LnSa.rock.Td1 <- f1.as(M, Rrup, T = Td1) + a12.as(T = Td1)*Frv +
      a13.as(T = Td1)*Fnm + a15.as(T = Td1)*Fas +
        f5.as(PGA1100,Vs30.rock, T = Td1) +
          Fhw*f4.as(Rjb, Rx, dip, Ztor, M, W, T = Td1) +
            f6.as(Ztor, T = Td1) + f8.as(Rrup, M, T = Td1) +
              f10.as(Z1.0.rock, Vs30.rock, T = Td1)
    # Upper bound Td for interpolation
    Td2 <- getPeriod(Td, "AS08")$upper
    LnSa.rock.Td2 <- f1.as(M, Rrup, T = Td2) + a12.as(T = Td2)*Frv +
      a13.as(T = Td2)*Fnm + a15.as(T = Td2)*Fas +
        f5.as(PGA1100, Vs30.rock, T = Td2) +
          Fhw*f4.as(Rjb, Rx, dip, Ztor, M, W, T = Td2) +
            f6.as(Ztor, T = Td2) + f8.as(Rrup, M, T = Td2) +
              f10.as(Z1.0.rock, Vs30.rock, T = Td2)
    # Actual Td, computed by interpolation
    LnSa.rock.Td <- interpolate(log(Td), log(Td1), log(Td2), LnSa.rock.Td1, LnSa.rock.Td2)

    # Calculate LnSa.rock.T (scale the Rock Sa to the spectral period T)
    LnSa.rock.T <- LnSa.rock.Td + log((Td/T)^2)
      
    # Calculate soil amplification
    SiteResponse.soil <- f5.as(PGA1100, Vs30, T) + f10.as(Z1.0, Vs30, T)
    SiteResponse.rock <- f5.as(PGA1100, Vs30.rock, T) +
                             f10.as(Z1.0.rock, Vs30.rock, T)
    Soil.amp <- SiteResponse.soil - SiteResponse.rock
      
    # Calculate final spectral acceleration
    LnSaConstDisp <- LnSa.rock.T + Soil.amp
     
    Sa <- exp(LnSaConstDisp)
  }
  return(Sa)
}




# 5. FINAL FUNCTION FOR AS08 GROUND MOTION CALCULATIONS
Sa.as <- function(M, Rjb, Vs30, VsFlag, epsilon, T, Rrup = NA, Rx = NA,
                  dip = NA, W = NA, Ztor = NA, Z1.0 = NA, rake = NA, Frv = NA,
                  Fnm = NA, Fhw = NA, azimuth = NA, Zhyp = NA, Fas = 0){

  # If T is a vector, perform calculation for each of the elements
  if(length(T) > 1) {
    return(sapply(T, Sa.as, M = M, Rjb = Rjb, Vs30 = Vs30, VsFlag = VsFlag,
                  epsilon = epsilon, Rrup = Rrup, Rx = Rx, dip = dip, W = W,
                  Ztor = Ztor, Z1.0 = Z1.0, rake = rake, Frv = Frv, Fnm = Fnm,
                  Fhw = Fhw, azimuth = azimuth, Zhyp = Zhyp, Fas = Fas))

    
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
    if(is.na(VsFlag == TRUE) | !(VsFlag == 1 | VsFlag == 0))
      stop("VsFlag must be either 1 (for measured Vs) or 0 (for estimated / inferred Vs)")
    if(is.na(epsilon) == TRUE)
      stop("epsilon must be numeric")
    if(is.na(Fas) == TRUE | !(Fas == 1 | Fas == 0))
      stop("Fas must be either 1 (for aftershocks) or 0 (for mainshocks, the default)")

    # Check optional distance, source, and site parameters
    if(is.na(Rrup) == FALSE & Rrup < 0)
      stop("Rrup must be a non-negative number; use NA if unknown.")
    if(is.na(dip) == FALSE & (dip <= 0 | dip > 90))
      stop("dip must be be greater than 0 and less than or equal to 90 deg; use NA if unknown.")
    if(is.na(W) == FALSE & W < 0)
      stop("W must be a non-negative number; use NA if unknown")
    if(is.na(Ztor) == FALSE & Ztor < 0)
      stop("Ztor must be a non-negative number; use NA if unknown.")
    if(is.na(Z1.0) == FALSE & Z1.0 < 0)
      stop("Z1.0 must be a non-negative number; use NA if unknown.")
    if(is.na(Zhyp) == FALSE & Zhyp < 0)
      stop("Zhyp must be a non-negative number; use NA if unknown.")

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
             "\n", "Fnm must be 1 for normal faulting (-120 <= rake <= -60) and 0 otherwise.")
    }
    # Ensure consistency between rake, Frv, and Fnm
    if(is.na(rake) == FALSE & is.na(Frv) == FALSE & is.na(Fnm) == FALSE){
      if(rake >= 30 & rake <= 150 & !(Frv == 1 & Fnm == 0))
        stop("Inconsistency between rake and style-of-faulting flag variables:", "\n",
             "Frv = 1 and Fnm = 0 for reverse faulting (30 <= rake <= 150)")
      if(rake >= -120 & rake <= -60 & !(Frv == 0 & Fnm == 1))
        stop("Inconsistency between rake and style-of-faulting flag variables:", "\n",
             "Frv = 0 and Fnm = 1 for normal faulting (-120 <= rake <= -60)")
      if(((rake < 30 & rake > -60) | rake < -120 | rake > 150) & !(Frv == 0 & Fnm == 0))
        stop("Inconsistency between rake and style-of-faulting flag variables:", "\n",
             "Frv = 0 and Fnm = 0 for strike-slip faulting.")
    }

   # Check hanging wall parameters
    if(is.na(azimuth) == TRUE & is.na(Fhw) == TRUE & is.na(Rx) == TRUE)
      stop("At least one of Rx, azimuth, and Fhw must be specified")
    if(is.na(azimuth) == FALSE & abs(azimuth) > 180)
      stop("the source-to-site azimuth must be between -180 and 180, inclusive")
    if(is.na(Fhw) == FALSE & !(Fhw == 1 | Fhw == 0))
      stop("Fhw must be either 1 (for sites on the hanging wall side of the fault)", "\n",
           "or 0 (for sites on the footwall side of the fault)")
    # Ensure consistency between Rx, azimuth, and Fhw
    if(is.na(dip) == FALSE & dip != 90){
      if(is.na(azimuth) == FALSE & is.na(Fhw) == FALSE){
        if(azimuth < 0 & Fhw == 1)
          stop("Inconsistency between azimuth and Fhw. Fhw must be 0 when azimuth < 0.")
        if(azimuth > 0 & Fhw == 0)
          stop("Inconsistency between azimuth and Fhw. Fhw must be 1 when azimuth > 0.")
      }
      if(is.na(Rx) == FALSE & is.na(Fhw) == FALSE){
        if(Rx < 0 & Fhw == 1)
          stop("Inconsistency between Rx and Fhw. Fhw must be 0 when Rx < 0.")
        if(Rx > 0 & Fhw == 0)
          stop("Inconsistency between Rx and Fhw. Fhw must be 1 when Rx > 0.")
      }
    }
    if(is.na(Rx) == FALSE & is.na(azimuth) == FALSE){
      if(!(Rx <= 0 & azimuth <= 0) & !(Rx >= 0 & azimuth >= 0))
        stop("Rx and azimuth must have the same sign.")
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
      if(rake >= -120 & rake <= -60)
        Fnm <- 1
      else
        Fnm <- 0
    }
    
    # Dip angle
    if(is.na(dip) == TRUE)
      dip <- dip.calc(rake)

    # Down-dip rupture width, W
    if(is.na(W) == TRUE)
      W <- W.calc(M, rake)
    
    # Depth to top of rupture, Ztor
    if(is.na(Ztor) == TRUE){
      # Hypocentral depth, Zhyp (used for calculating Ztor)
      if(is.na(Zhyp) == TRUE)
        Zhyp <- Zhyp.calc(M, rake)
      # Calculation of Ztor
      Ztor <- Ztor.calc(W, dip, Zhyp)
    }

    # Determine hanging wall flag, Fhw
      # Calculate from azimuth if provided
      if(is.na(azimuth) == FALSE){
        if(azimuth > 0 & azimuth < 180 & dip != 90)
          Fhw <- 1
        else
          Fhw <- 0
      # Calculate from Rx if azimuth is not provided
      } else{
        if(is.na(Rx) == FALSE){
          if(Rx > 0 & dip != 90)
            Fhw <- 1
          else
            Fhw <- 0
        }
      }
      # For the special case of vertical fault, ensure that Fhw = 0
      if(dip == 90)
        Fhw <- 0

    # Azimuth angle
    if(is.na(azimuth) == TRUE){
      # If hanging wall site, assume azimuth = 50 deg
      if(Fhw == 1)
        azimuth <- 50
      # If footwall site, assume azimuth = -50 deg  
      else if(Fhw == 0)
        azimuth <- -50
    }
    # Ensure that azimuth = 90 for Rjb = 0
    if(Rjb == 0){
      Fhw <- 1
      azimuth <- 90
    }
      
    # Site coordinate, Rx
    if(is.na(Rx) == TRUE)
      Rx <- Rx.calc(Rjb, Ztor, W, dip, azimuth, Rrup)
  
    # Rupture distance, Rrup
    if(is.na(Rrup) == TRUE)
      Rrup <- Rrup.calc(Rx, Ztor, W, dip, azimuth, Rjb)

    # Depth parameter, Z1.0
    if(is.na(Z1.0) == TRUE)
      Z1.0 <- Z1.calc.as(Vs30)

    
    # 5C. CALCULATE GROUND MOTION PARAMETER
    
    # Is interpolation necessary?
    interp <- getPeriod(T, "AS08")$interp
    
    # If interpolation is not necessary, compute Sa
    if(interp == FALSE){
      LnSaMedian <- log(SaMedian.as(M, Rjb, Rrup, Rx, Ztor, Frv, Fnm,
                                    Fas, Fhw, dip, Vs30, Z1.0, W, T))
      PGA1100 <- PGA1100.as(M, Rrup, Rjb, Rx, Ztor, Frv, Fnm, Fas, Fhw, dip, W)  
      epsilon.sigmaTot <- epsilon * SigmaTot.as(M, PGA1100, Vs30, VsFlag, T)
      LnSa <- LnSaMedian + epsilon.sigmaTot
      return(exp(LnSa))
    } else{

    # If interpolation is necessary, compute Sa
      if(interp == TRUE){
        T1 <- getPeriod(T, "AS08")$lower
        T2 <- getPeriod(T, "AS08")$upper
        
        # PGA1100
        PGA1100 <- PGA1100.as(M, Rrup, Rjb, Rx, Ztor, Frv, Fnm, Fas, Fhw, dip, W)    
      
        # Calculation for T1
        LnSaMedian.T1 <- log(SaMedian.as(M, Rjb, Rrup, Rx, Ztor, Frv, Fnm,
                                         Fas, Fhw, dip, Vs30, Z1.0, W, T1))
        epsilon.sigmaTot.T1 <- epsilon * SigmaTot.as(M, PGA1100, Vs30, VsFlag, T1)
        LnSaT1 <- LnSaMedian.T1 + epsilon.sigmaTot.T1
      
        # Calculation for T2
        LnSaMedian.T2 <- log(SaMedian.as(M, Rjb, Rrup, Rx, Ztor, Frv, Fnm,
                                         Fas, Fhw, dip, Vs30, Z1.0, W, T2))
        epsilon.sigmaTot.T2 <- epsilon * SigmaTot.as(M, PGA1100, Vs30, VsFlag, T2)
        LnSaT2 <- LnSaMedian.T2 + epsilon.sigmaTot.T2
      
        # Interpolated value
        LnSa <- interpolate(log(T), log(T1), log(T2), LnSaT1, LnSaT2)
        return(exp(LnSa))
      }
    }
  }
}

