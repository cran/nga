# BOORE & ATKINSON NGA MODEL
# Boore, D. M., and G. M. Atkinson (2008). Ground-Motion Prediction Equations for Average
# Horizontal Component of PGA, PGV, and PSA. Earthquake Spectra, Vol. 24, pp. 99-138.


# OUTLINE OF CODE
# 1. Model Coefficients
#    a. Periods with defined coefficients
#    b. Period-independent coefficients
#    c. Period-dependent coefficients for median ground motion term
#    d. Coefficients for standard deviation term
# 2. Necessary Functions for Calculating Median Ground Motion
#    a. Distance Function
#    b. Magnitude Function
#    c. Site Amplification Function
# 3. Median Ground Motion Calculation of Sa, PGA, and PGV
# 4. Final Function for BA08 Ground Motion Calculations
#    a. Check input parameters
#    b. Obtain estimates of unspecified input parameters
#    c. Calculate ground motion parameter
#    **** This is the function that users will primarily use for ground motion calculations ****


# INPUT PARAMETERS FOR THE FUNCTIONS BELOW:
#   M = Moment magnitude
#   Rjb = Joyner-Boore distance (km)
#   Vs30 = Time-averaged shear wave velocity over 30 m subsurface depth  (m/sec)
#   rake = Rake angle of fault movement (deg)
#   epsilon = number of standard deviations to be considered in the calculations
#   T = Spectral period, sec (0 for PGA; -1 for PGV)

# OUTPUT PARAMETERS (from Sa function):
#   Sa = Spectral acceleration (g)
#   PGA = Peak ground acceleration (g); calculated by evaluating Sa at T = 0
#   PGV = Peak ground velocity (cm/sec); calculated by evaluating Sa at T = -1




# 1. MODEL COEFFICIENTS

# 1a. Periods with defined coefficients (PGA is 0; PGV is -1)
periods.ba <- function(positive = FALSE){

  # Return list of periods excluding PGA and PGV
  if(positive){
    T <- c(0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40,
           0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0)

  # Return list of periods including PGA and PGV
  }else{
    T <- c(0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40,
           0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0, 0, -1)
  }
  return(T)
}


# 1b. Period-independent coefficients
coefs.ba <- function(){
  
  # Site-amplification coefficients
  # Table 4 in Boore and Atkinson (2008)
  c("a1 = 0.03", "pga.low = 0.06", "a2 = 0.09", "V1 = 180", "V2 = 300", "Vref = 760",

  # Distance-scaling coefficients
  # Table 6 in Boore and Atkinson (2008)
    "Mref = 4.5", "Rref = 1.0")
}



# 1c. Period-dependent coefficients for median ground motion term

# Site amplification coefficients
# Table 3 in Boore and Atkinson (2008)

blin.ba <- function(T) {
  Period.list <- periods.ba()
  blin.list <- c(-0.36, -0.34, -0.33, -0.29, -0.23, -0.25, -0.28, -0.31, -0.39,
                 -0.44, -0.5, -0.6, -0.69, -0.7, -0.72, -0.73, -0.74, -0.75,
                 -0.75, -0.692, -0.65, -0.36, -0.6)
  blin.list[match(T, Period.list)]
}

b1.ba <- function(T) {
  Period.list <- periods.ba()
  b1.list <- c(-0.64, -0.63, -0.62, -0.64, -0.64, -0.6, -0.53, -0.52, -0.52,
               -0.52, -0.51, -0.5, -0.47, -0.44, -0.4, -0.38, -0.34, -0.31,
               -0.291, -0.247, -0.215, -0.64, -0.5)
  b1.list[match(T, Period.list)]
}

b2.ba <- function(T) {
  Period.list <- periods.ba()
  b2.list <- c(-0.14, -0.12, -0.11, -0.11, -0.11, -0.13, -0.18, -0.19, -0.16,
               -0.14, -0.1, -0.06, 0 , 0, 0, 0, 0, 0, 0, 0, 0, -0.14, -0.06)
  b2.list[match(T, Period.list)]
}


# Distance-scaling coefficients
# Table 6 in Boore and Atkinson (2008)

c1.ba <- function(T) {
  Period.list <- periods.ba()
  c1.list <- c(-0.66220, -0.66600, -0.69010, -0.71700, -0.72050, -0.70810,
               -0.69610, -0.58300, -0.57260, -0.55430, -0.64430, -0.69140,
               -0.74080, -0.81830, -0.83030, -0.82850, -0.78440, -0.68540,
               -0.50960, -0.37240, -0.09824, -0.66050, -0.87370)
  c1.list[match(T, Period.list)]
}

c2.ba <- function(T) {
  Period.list <- periods.ba()
  c2.list <- c(0.12000, 0.12280, 0.12830, 0.13170, 0.12370, 0.11170, 0.09884,
               0.04273, 0.02977, 0.01955, 0.04394, 0.06080, 0.07518, 0.10270,
               0.09793, 0.09432, 0.07282, 0.03758, -0.02391, -0.06568, -0.13800,
               0.11970, 0.10060)  
  c2.list[match(T, Period.list)]
}

c3.ba <- function(T) {
  Period.list <- periods.ba()
  c3.list <- c(-0.01151, -0.01151, -0.01151, -0.01151, -0.01151, -0.01151,
               -0.01113, -0.00952, -0.00837, -0.00750, -0.00626, -0.00540,
               -0.00409, -0.00334, -0.00255, -0.00217, -0.00191, -0.00191,
               -0.00191, -0.00191, -0.00191, -0.01151, -0.00334)  
  c3.list[match(T, Period.list)]
}

h.ba <- function(T) {
  Period.list <- periods.ba()
  h.list <- c(1.35, 1.35, 1.35, 1.35, 1.55, 1.68, 1.86, 1.98, 2.07, 2.14, 2.24, 2.32,
              2.46, 2.54, 2.66, 2.73, 2.83, 2.89, 2.93, 3.00, 3.04, 1.35, 2.54)
  h.list[match(T, Period.list)]
}


# Magnitude-scaling coefficients
# Table 7 in Boore and Atkinson (2008)

e1.ba <- function(T) {
  Period.list <- periods.ba()
  e1.list <- c(-0.52883, -0.52192, -0.45285, -0.28476, 0.00767, 0.20109, 0.46128,
               0.57180, 0.51884, 0.43825, 0.39220, 0.18957, -0.21338, -0.46896,
               -0.86271, -1.22652, -1.82979, -2.24656, -1.28408, -1.43145,
               -2.15446, -0.53804, 5.00121)
  e1.list[match(T, Period.list)]
}

e2.ba <- function(T) {
  Period.list <- periods.ba()
  e2.list <- c(-0.49429, -0.48508, -0.41831, -0.25022, 0.04912, 0.23102, 0.48661,
               0.59253, 0.53496, 0.44516, 0.40602, 0.19878, -0.19496, -0.43443,
               -0.79593, -1.15514, -1.74690, -2.15906, -1.21270, -1.31632,
               -2.16137, -0.50350, 5.04727)
  e2.list[match(T, Period.list)]
}

e3.ba <- function(T) {
  Period.list <- periods.ba()
  e3.list <- c(-0.74551, -0.73906, -0.66722, -0.48462, -0.20578, 0.03058, 0.30185,
               0.4086, 0.3388, 0.25356, 0.21398, 0.00967, -0.49176, -0.78465,
               -1.20902, -1.57697, -2.22584, -2.58228, -1.50904, -1.81022,
               -2.53323, -0.75472, 4.63188)
  e3.list[match(T, Period.list)]
}

e4.ba <- function(T) {
  Period.list <- periods.ba()
  e4.list <- c(-0.49966, -0.48895, -0.42229, -0.26092, 0.02706, 0.22193, 0.49328,
               0.61472, 0.57747, 0.5199, 0.4608, 0.26337, -0.10813, -0.3933,
               -0.88085, -1.27669, -1.91814, -2.38168, -1.41093, -1.59217,
               -2.14635, -0.5097, 5.0821)
  e4.list[match(T, Period.list)]
}

e5.ba <- function(T) {
  Period.list <- periods.ba()
  e5.list <- c(0.28897, 0.25144, 0.17976, 0.06369, 0.0117, 0.04697, 0.1799, 0.52729,
               0.6088, 0.64472, 0.7861, 0.76837, 0.75179, 0.6788, 0.70689, 0.77989,
               0.77966, 1.24961, 0.14271, 0.52407, 0.40387, 0.28805, 0.18322)
  e5.list[match(T, Period.list)]
}

e6.ba <- function(T) {
  Period.list <- periods.ba()
  e6.list <- c(-0.10019, -0.11006, -0.12858, -0.15752, -0.17051, -0.15948,
               -0.14539, -0.12964, -0.13843, -0.15694, -0.07843, -0.09054,
               -0.14053, -0.18257, -0.2595, -0.29657, -0.45384, -0.35874,
               -0.39006, -0.37578, -0.48492, -0.10164, -0.12736)
  e6.list[match(T, Period.list)]
}

e7.ba <- function(T) {
  Period.list <- periods.ba()
  e7.list <- c(0, 0, 0, 0, 0, 0, 0, 0.00102, 0.08607, 0.10601, 0.02262, 0, 0.10302,
               0.05393, 0.19082, 0.29888, 0.67466, 0.79508, 0, 0, 0, 0, 0)
  e7.list[match(T, Period.list)]
}

Mh.ba <- function(T) {
  Period.list <- periods.ba()
  Mh.list <- c(6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75,
               6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75,
               8.50, 8.50, 8.50, 6.75, 8.50)
  Mh.list[match(T, Period.list)]
}




# 1d. Coefficients for standard deviation term
#     Table 8 in Boore and Atkinson (2008)

# Intra-event SD
Sigma.ba <- function(T) {
  Period.list <- periods.ba()
  Sigma.list <- c(0.502, 0.502, 0.507, 0.516, 0.513, 0.520, 0.518, 0.523, 0.527,
                  0.546, 0.541, 0.555, 0.571, 0.573, 0.566, 0.580, 0.566, 0.583,
                  0.601, 0.626, 0.645, 0.502, 0.500)
  Sigma.list[match(T, Period.list)]
}

# Inter-event SD when the fault mechanism is unspecified
TauU.ba <- function(T) { 
  Period.list <- periods.ba()
  TauU.list <- c(0.267, 0.267, 0.276, 0.286, 0.322, 0.313, 0.288, 0.283, 0.267,
                 0.272, 0.267, 0.265, 0.311, 0.318, 0.382, 0.398, 0.410, 0.394,
                 0.414, 0.465, 0.355, 0.265, 0.286)
  TauU.list[match(T, Period.list)]
}

# Inter-event SD when the fault mechanism is specified
TauM.ba <- function(T) { 
  Period.list <- periods.ba()
  TauM.list <- c(0.262, 0.262, 0.274, 0.286, 0.320, 0.318, 0.290, 0.288, 0.267,
                 0.269, 0.267, 0.265, 0.299, 0.302, 0.373, 0.389, 0.401, 0.385,
                 0.437, 0.477, 0.477, 0.260, 0.256)
  TauM.list[match(T, Period.list)]
}

# Total SD when the fault mechanism is unspecified
# Calculated by sqrt(Sigma^2 + TauU^2)
SigmaTotU.ba <- function(T) { 
  Period.list <- periods.ba()
  SigmaTotU.list <-  c(0.569, 0.569, 0.578, 0.589, 0.606, 0.608, 0.592, 0.596, 0.592,
                       0.608, 0.603, 0.615, 0.649, 0.654, 0.684, 0.702, 0.700, 0.702,
                       0.730, 0.781, 0.735, 0.566, 0.576)
  SigmaTotU.list[match(T, Period.list)]
}

# Total SD when the fault mechanism is specified
# Calculated by sqrt(Sigma^2 + TauM^2)
SigmaTotM.ba <- function(T) { 
  Period.list <- periods.ba()
  SigmaTotM.list <- c(0.566, 0.566, 0.576, 0.589, 0.606, 0.608, 0.594, 0.596, 0.592,
                      0.608, 0.603, 0.615, 0.645, 0.647, 0.679, 0.700, 0.695, 0.698,
                      0.744, 0.787, 0.801, 0.564, 0.560)
  SigmaTotM.list[match(T, Period.list)]
}




# 2. NECESSARY FUNCTIONS FOR CALCULATING MEDIAN GROUND MOTION

# 2a. Distance Function
Fd.ba <- function(M, Rjb, T)
  { 
    # Load period-independent coefficients
    coefs.list <- coefs.ba()
    for(i in 1:length(coefs.list))
      eval(parse(text = coefs.list[i]))
 
    # Calculate R (Eqn 3)
    R <- sqrt(Rjb^2 + (h.ba(T))^2)

    # Return Fd (Eqn 4)
    return((c1.ba(T) + c2.ba(T)*(M - Mref))*log(R / Rref) + c3.ba(T)*(R - Rref))
  }


# 2b. Magnitude Function
Fm.ba <- function(M, U, SS, NS, RS, T)
  {
    if(M <= Mh.ba(T))  # Eqn 5A
      {
        return(e1.ba(T)*U + e2.ba(T)*SS + e3.ba(T)*NS +
               e4.ba(T)*RS + e5.ba(T)*(M - Mh.ba(T)) + e6.ba(T)*(M - Mh.ba(T))^2)
      } else
    if(M > Mh.ba(T))  # Eqn 5B
      {
        return(e1.ba(T)*U + e2.ba(T)*SS + e3.ba(T)*NS +
               e4.ba(T)*RS + e7.ba(T)*(M - Mh.ba(T)))
      }
  }


# 2c. Site Amplification Function
Fs.ba <- function(M, Rjb, Vs30, U, SS, NS, RS, T)
  {

    # Load period-independent coefficients
    coefs.list <- coefs.ba()
    for(i in 1:length(coefs.list))
      eval(parse(text = coefs.list[i]))

    
    # i. LINEAR TERM (Eqn 7)
    Flin <- blin.ba(T)*log(Vs30 / Vref)

    
    # ii. NONLINEAR TERM
    
    # Calculate pga4nl, the median PGA when Vs30 = Vref = 760 m/s
    pga4nl <- exp(Fm.ba(M, U, SS, NS, RS, T=0) + Fd.ba(M, Rjb, T=0))

    # Calculate Nonlinear slope, bnl (Eqns 13a, 13b, and 13c)
    if(Vs30 <= V1) {
      bnl <- b1.ba(T)
    } else {
      if(Vs30 > V1 & Vs30 <= V2) {
        bnl <- (b1.ba(T) - b2.ba(T))*log(Vs30/V2) / log(V1/V2) + b2.ba(T)
      } else { 
        if(Vs30 > V2 & Vs30 < Vref) {
          bnl <- b2.ba(T)*log(Vs30/Vref) / log(V2/Vref)
        } else {
          if(Vref <= Vs30) {
            bnl <- 0
          }
        }
      }
    }

    # Calculate smoothing constants (Eqns 9, 10, 11, and 12)
    dx <- log(a2/a1)
    dy <- bnl*log(a2/pga.low)
    c <- (3*dy - bnl*dx)/(dx^2)
    d <- -(2*dy - bnl*dx)/(dx^3)
  
    # Final equation for nonlinear term (Eqns 8a, 8b, and 8c)
    if(pga4nl <= a1) {
      Fnl <- bnl*log(pga.low/0.1)
    } else {
      if(pga4nl > a1 & pga4nl <= a2) {
        Fnl <- bnl*log(pga.low/0.1) + c*(log(pga4nl/a1))^2 + d*(log(pga4nl/a1))^3
      } else {
        if(pga4nl > a2) {
          Fnl <- bnl*log(pga4nl/0.1)
        }
      }
    }

    # Return Fs (Eqn 6)
    return(Flin + Fnl)
  }




# 3. MEDIAN GROUND MOTION CALCULATION of Sa, PGA, and PGV  (Eqn 1)
SaMedian.ba <- function(M, Rjb, Vs30, U, SS, NS, RS, T) {
  SaMedian <- exp(Fm.ba(M, U, SS, NS, RS, T) + Fd.ba(M, Rjb, T) +
                  Fs.ba(M, Rjb, Vs30, U, SS, NS, RS, T))
  return(SaMedian)
}




# 4. FINAL FUNCTION FOR BA08 GROUND MOTION CALCULATIONS
Sa.ba <- function(M, Rjb, Vs30, rake = NA, U = NA, SS = NA, NS = NA, RS = NA, epsilon, T){

  # If T is a vector, perform calculation for each of the elements
  if(length(T) > 1) {
    return(sapply(T, Sa.ba, M = M, Rjb = Rjb, Vs30 = Vs30, rake = rake,
                  U = U, SS = SS, NS = NS, RS = RS, epsilon = epsilon))

  # Perform calculation for single value of T:
  } else {

    
    # 4A. CHECK INPUT PARAMETERS

    # Check mandatory input parameters
    if(is.na(M) == TRUE | M < 0)
      stop("M must be a positive number")
    if(is.na(Rjb) == TRUE | Rjb < 0)
      stop("Rjb must be a non-negative number")
    if(is.na(Vs30) == TRUE | Vs30 < 0)
      stop("Vs30 must be a positive number")
    if(is.na(epsilon) == TRUE)
      stop("epsilon must be numeric")

    # Check style-of-faulting parameters
    if(is.na(rake) == TRUE & (is.na(U) == TRUE | is.na(SS) == TRUE |
                   is.na(NS) == TRUE | is.na(RS) == TRUE))
      stop("either (1) the rake angle, or (2) all the style-of-faulting flag variables", "\n",
           "(U, SS, NS, and RS) must be specified")
    if(is.na(rake) == FALSE & (is.na(U) == FALSE | is.na(SS) == FALSE |
                   is.na(NS) == FALSE | is.na(RS) == FALSE) &
       (is.na(U) == TRUE | is.na(SS) == TRUE | is.na(NS) == TRUE |
        is.na(RS) == TRUE))
      stop("either (1) the rake angle, or (2) all the style-of-faulting flag variables", "\n",
           "(U, SS, NS, and RS) must be specified")
    if(is.na(rake) == TRUE & U + SS + NS + RS != 1)
      stop("exactly one style-of-faulting flag variable (U, SS, NS, or RS) must be equal to 1")
  
    # Check for consistency, if both rake and style-of-faulting parameters are supplied
    if(is.na(rake) == FALSE & (is.na(U) == FALSE | is.na(SS) == FALSE |
                   is.na(NS) == FALSE | is.na(RS) == FALSE)){
      if(U + SS + NS + RS != 1)
        stop("exactly one style-of-faulting flag variable (U, SS, NS, or RS) must be equal to 1")
      if(U == 1)
        stop("for unspecified faulting (U = 1), the rake should not be passed as an argument")
      else{
        if(NS == 1 & !(rake >= -150 & rake <= -30))
          stop("inconsistency between rake and NS:  for NS = 1, the rake must be between -150 and -30, inclusive")
        if(RS == 1 & !(rake >=30 & rake <= 150))
          stop("inconsistency between rake and RS:  for RS = 1, the rake must be between 30 and 150, inclusive")
        if(SS == 1 & !(abs(rake) < 30 | abs(rake) > 150))
          stop("inconsistency between rake and SS:  for SS = 1, the rake cannot be in the ranges for normal or reverse faulting")
      }
    }


    # 4B. OBTAIN ESTIMATES OF UNSPECIFIED INPUT PARAMETERS
    
    # Convert rake to fault type, if fault type not provided
    if(is.na(U) == TRUE){
      # Normal faulting
      if(rake >= -150 & rake <= -30){
        U <- 0
        NS <- 1
        RS <- 0
        SS <- 0
        # Reverse faulting   
      } else{
        if(rake >= 30 & rake <= 150){
          U <- 0
          NS <- 0
          RS <- 1
          SS <- 0
          # Strike-slip faulting
        } else{
          U <- 0
          NS <- 0
          RS <- 0
          SS <- 1
        }
      }
    }


    # 4C. CALCULATE GROUND MOTION PARAMETER
    
    # Is interpolation necessary?
    interp <- getPeriod(T, "BA08")$interp

    # If interpolation is not necessary, compute Sa
    if(interp == FALSE){
      # Compute median
      LnSaMedian <- log(SaMedian.ba(M, Rjb, Vs30, U, SS, NS, RS, T))
      # Compute SD, depending on whether or not the fault mechanism is specified
      if(U == 0)
        epsilon.sigmaTot <- epsilon * SigmaTotM.ba(T)   # Specified fault mechanism
      else
        epsilon.sigmaTot <- epsilon * SigmaTotU.ba(T)   # Unspecified fault mechanism
      LnSa <- LnSaMedian + epsilon.sigmaTot
      return(exp(LnSa))
    } else{


    # If interpolation is necessary, compute Sa
      if(interp == TRUE){
        T1 <- getPeriod(T, "BA08")$lower
        T2 <- getPeriod(T, "BA08")$upper    
      
        # Calculation for T1
        LnSaMedian.T1 <- log(SaMedian.ba(M, Rjb, Vs30, U, SS, NS, RS, T1))
        if(U == 0)
          epsilon.sigmaTot.T1 <- epsilon * SigmaTotM.ba(T1)   # Specified fault mechanism
        else
          epsilon.sigmaTot.T1 <- epsilon * SigmaTotU.ba(T1)   # Unspecified fault mechanism
        LnSaT1 <- LnSaMedian.T1 + epsilon.sigmaTot.T1
      
        # Calculation for T2
        LnSaMedian.T2 <- log(SaMedian.ba(M, Rjb, Vs30, U, SS, NS, RS, T2))
        if(U == 0)
          epsilon.sigmaTot.T2 <- epsilon * SigmaTotM.ba(T2)   # Specified fault mechanism
        else
          epsilon.sigmaTot.T2 <- epsilon * SigmaTotU.ba(T2)   # Unspecified fault mechanism
        LnSaT2 <- LnSaMedian.T2 + epsilon.sigmaTot.T2
      
        # Interpolated value
        LnSa <- interpolate(log(T), log(T1), log(T2), LnSaT1, LnSaT2)
        return(exp(LnSa))
      }
    }
  }
}

