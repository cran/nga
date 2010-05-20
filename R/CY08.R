# CHIOU & YOUNGS NGA MODEL
# Chiou, B. S.-J., and R. R. Youngs (2008). An NGA Model for the Average
# Horizontal Component of Peak Ground Motion and Response Spectra.
# Earthquake Spectra, Vol. 24, pp. 67-97.


# OUTLINE OF CODE
# 1. Model Coefficients
#    a. Periods with defined coefficients
#    b. Period-independent coefficients
#    c. Period-dependent coefficients for median ground motion term
#    d. Coefficients for standard deviation term
# 2. Necessary Function for Calculating Median Ground Motion
#    Ground motion on reference rock site
# 3. Necessary Functions for Calculating Standard Deviation Term
#    a. Nonlinear site response term
#    b. Sigma, intra-event standard deviation
#    c. Tau, inter-event standard deviation
#    d. SigmaTot, total standard deviation
# 4. Median Ground Motion Calculation of Sa, PGA, and PGV
#    Site response function
# 5. Final Function for CY08 Ground Motion Calculations
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
#   Fhw = Hanging wall flag (1 for site on hanging wall site of fault,
#         0 otherwise); calculated from azimuth or Rx (if provided)
#   azimuth = source-to-site azimuth (deg); see Figure 2 in Kaklamanos and Baise (2010)
#   Zhyp = hypocentral depth (km)
#   AS = Aftershock flag (1 for aftershocks, 0 for mainshocks)
#   Yref = median ground motion for the reference rock site condition (Eqn 13a)


# OUTPUT PARAMETERS (from Sa function):
#   Sa = Spectral acceleration (g)
#   PGA = Peak ground acceleration (g); calculated by evaluating Sa at T = 0
#   PGV = Peak ground velocity (cm/sec); calculated by evaluating Sa at T = -1




# 1. MODEL COEFFICIENTS

# 1a. Periods with defined coefficients (PGA is 0; PGV is -1)
periods.cy <- function(positive = FALSE) {
  
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
#     Table 1 in Chiou and Youngs (2008)
coefs.cy <- function() {
  c("c2 = 1.06", "c3 = 3.45", "c4 = -2.1", "c4a = -0.5", "cRB = 50", "cHM = 3", "cGamma3 = 4")
}


# 1c. Period-dependent coefficients for median ground motion term
#     Tables 2 and 3 in Chiou and Youngs (2008)

# Table 2 - Coefficients of Model for Ln(Yref)

c1.cy <- function(T) {
  Period.list <- periods.cy()
  c1.list <- c(-1.2687, -1.2515, -1.1744, -1.0671, -0.9464, -0.7051, -0.5747,
               -0.5309, -0.6352, -0.7766, -0.9278, -1.2176, -1.4695, -1.9278,
               -2.2453, -2.7307, -3.1413, -3.7413, -4.1814, -4.5187, -5.1224,
               -5.5872, -1.2687, 2.2884)
  c1.list[match(T, Period.list)]
}

c1a.cy <- function(T) {
  Period.list <- periods.cy()
  c1a.list <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.0999,
                0.0997, 0.0991, 0.0936, 0.0766, 0.0022, -0.0591, -0.0931, -0.0982,
                -0.0994, -0.0999, -0.1, 0.1, 0.1094)
  c1a.list[match(T, Period.list)]
}

c1b.cy <- function(T) {
  Period.list <- periods.cy()
  c1b.list <- c(-0.255, -0.255, -0.255, -0.255, -0.255, -0.254, -0.253, -0.25,
                -0.2449, -0.2382, -0.2313, -0.2146, -0.1972, -0.162, -0.14, -0.1184,
                -0.11, -0.104, -0.102, -0.101, -0.101, -0.1, -0.255, -0.0626)
  c1b.list[match(T, Period.list)]
}

cn.cy <- function(T) {
  Period.list <- periods.cy()
  cn.list <- c(2.996, 3.292, 3.514, 3.563, 3.547, 3.448, 3.312, 3.044, 2.831,
               2.658, 2.505, 2.261, 2.087, 1.812, 1.648, 1.511, 1.47, 1.456,
               1.465, 1.478, 1.498, 1.502, 2.996, 1.648)
  cn.list[match(T, Period.list)]
}

cM.cy <- function(T) {
  Period.list <- periods.cy()
  cM.list <- c(4.1840, 4.1879, 4.1556, 4.1226, 4.1011, 4.0860, 4.1030, 4.1717,
               4.2476, 4.3184, 4.3844, 4.4979, 4.5881, 4.7571, 4.8820, 5.0697,
               5.2173, 5.4385, 5.5977, 5.7276, 5.9891, 6.1930, 4.1840, 4.2979)
  cM.list[match(T, Period.list)]
}

c5.cy <- function(T) {
  Period.list <- periods.cy()
  c5.list <- c(6.16, 6.158, 6.155, 6.1508, 6.1441, 6.12, 6.085, 5.9871, 5.8699,
               5.7547, 5.6527, 5.4997, 5.4029, 5.29, 5.248, 5.2194, 5.2099,
               5.204, 5.202, 5.201, 5.2, 5.2, 6.16, 5.17)
  c5.list[match(T, Period.list)]
}

c6.cy <- function(T) {
  Period.list <- periods.cy()
  c6.list <- c(0.4893, 0.4892, 0.489, 0.4888, 0.4884, 0.4872, 0.4854, 0.4808,
               0.4755, 0.4706, 0.4665, 0.4607, 0.4571, 0.4531, 0.4517, 0.4507,
               0.4504, 0.4501, 0.4501, 0.45, 0.45, 0.45, 0.4893, 0.4407)
  c6.list[match(T, Period.list)]
}

c7.cy <- function(T) {
  Period.list <- periods.cy()
  c7.list <- c(0.0512, 0.0512, 0.0511, 0.0508, 0.0504, 0.0495, 0.0489, 0.0479,
               0.0471, 0.0464, 0.0458, 0.0445, 0.0429, 0.0387, 0.035, 0.028,
               0.0213, 0.0106, 0.0041, 0.001, 0, 0, 0.0512, 0.0207)
  c7.list[match(T, Period.list)]
}

c7a.cy <- function(T) {
  Period.list <- periods.cy()
  c7a.list <- c(0.0860, 0.0860, 0.0860, 0.0860, 0.0860, 0.0860, 0.0860,
                0.0860, 0.0860, 0.0860, 0.0860, 0.0850, 0.0830, 0.0690,
                0.0450, 0.0134, 0.0040, 0.0010, 0, 0, 0, 0, 0.0860, 0.0437)
  c7a.list[match(T, Period.list)]
}

c9.cy <- function(T) {
  Period.list <- periods.cy()
  c9.list <- c(0.79, 0.8129, 0.8439, 0.874, 0.8996, 0.9442, 0.9677, 0.966,
               0.9334, 0.8946, 0.859, 0.8019, 0.7578, 0.6788, 0.6196,
               0.5101, 0.3917, 0.1244, 0.0086, 0, 0, 0, 0.79, 0.3079)
  c9.list[match(T, Period.list)]
}

c9a.cy <- function(T) {
  Period.list <- periods.cy()
  c9a.list <- c(1.5005, 1.5028, 1.5071, 1.5138, 1.523, 1.5597, 1.6104, 1.7549,
                1.9157, 2.0709, 2.2005, 2.3886, 2.5, 2.6224, 2.669, 2.6985,
                2.7085, 2.7145, 2.7164, 2.7172, 2.7177, 2.718, 1.5005, 2.669)
  c9a.list[match(T, Period.list)]
}

c10.cy <- function(T) {
  Period.list <- periods.cy()
  c10.list <-  c(-0.3218, -0.3323, -0.3394, -0.3453, -0.3502, -0.3579,
                 -0.3604, -0.3565, -0.3470, -0.3379, -0.3314, -0.3256,
                 -0.3189, -0.2702, -0.2059, -0.0852, 0.0160, 0.1876,
                 0.3378, 0.4579, 0.7514, 1.1856, -0.3218, -0.1166)
  c10.list[match(T, Period.list)]
}

cGamma1.cy <- function(T) {
  Period.list <- periods.cy()
  cGamma1.list <- c(-0.00804, -0.00811, -0.00839, -0.00875, -0.00912, -0.00973,
                     -0.00975, -0.00883, -0.00778, -0.00688, -0.00612, -0.00498,
                     -0.00420, -0.00308, -0.00246, -0.00180, -0.00147, -0.00117,
                     -0.00107, -0.00102, -0.00096, -0.00094, -0.00804, -0.00275)
  cGamma1.list[match(T, Period.list)]
}

cGamma2.cy <- function(T) {
  Period.list <- periods.cy()
  cGamma2.list <- c(-0.00785, -0.00792, -0.00819, -0.00855, -0.00891, -0.00950,
                     -0.00952, -0.00862, -0.00759, -0.00671, -0.00598, -0.00486,
                     -0.00410, -0.00301, -0.00241, -0.00176, -0.00143, -0.00115,
                     -0.00104, -0.00099, -0.00094, -0.00091, -0.00785, -0.00625)
  cGamma2.list[match(T, Period.list)]
}


# Table 3 - Coefficients of site response model for Ln(Y)

phi1.cy <- function(T) {
  Period.list <- periods.cy()
  phi1.list <- c(-0.4417, -0.434, -0.4177, -0.4, -0.3903, -0.404, -0.4423,
                 -0.5162, -0.5697, -0.6109, -0.6444, -0.6931, -0.7246, -0.7708,
                 -0.799, -0.8382, -0.8663,  -0.9032, -0.9231, -0.9222, -0.8346,
                 -0.7332, -0.4417, -0.7861)
  phi1.list[match(T, Period.list)]
}

phi2.cy <- function(T) {
  Period.list <- periods.cy()
  phi2.list <- c(-0.1417, -0.1364, -0.1403, -0.1591, -0.1862, -0.2538, -0.2943,
                 -0.3113, -0.2927, -0.2662, -0.2405, -0.1975, -0.1633, -0.1028,
                 -0.0699, -0.0425, -0.0302, -0.0129, -0.0016, 0, 0, 0, -0.1417,
                 -0.0699)
  phi2.list[match(T, Period.list)]
}

phi3.cy <- function(T) {
  Period.list <- periods.cy()
  phi3.list <- c(-0.00701, -0.007279, -0.007354, -0.006977, -0.006467,
                 -0.005734, -0.005604, -0.005845, -0.006141, -0.006439,
                 -0.006704, -0.007125, -0.007435, -0.00812, -0.008444,
                 -0.007707, -0.004792, -0.001828, -0.001523, -0.00144,
                 -0.001369, -0.001361, -0.00701, -0.008444)
  phi3.list[match(T, Period.list)]
}

phi4.cy <- function(T) {
  Period.list <- periods.cy()
  phi4.list <- c(0.102151, 0.10836, 0.119888, 0.133641, 0.148927, 0.190596,
                 0.230662, 0.266468, 0.255253, 0.231541, 0.207277, 0.165464,
                 0.133828, 0.085153, 0.058595, 0.031787, 0.019716, 0.009643,
                 0.005379, 0.003223, 0.001134, 0.000515, 0.102151, 5.41)
  phi4.list[match(T, Period.list)]
}

phi5.cy <- function(T) {
  Period.list <- periods.cy()
  phi5.list <- c(0.2289, 0.2289, 0.2289, 0.2289, 0.229, 0.2292, 0.2297, 0.2326,
                 0.2386, 0.2497, 0.2674, 0.312, 0.361, 0.4353, 0.4629, 0.4756,
                 0.4785, 0.4796, 0.4799, 0.4799, 0.48, 0.48, 0.2289, 0.2899)
  phi5.list[match(T, Period.list)]
}

phi6.cy <- function(T) {
  Period.list <- periods.cy()
  phi6.list <- c(0.014996, 0.014996, 0.014996, 0.014996, 0.014996, 0.014996,
                 0.014996, 0.014988, 0.014964, 0.014881, 0.014639, 0.013493,
                 0.011133, 0.006739, 0.005749, 0.005544, 0.005521, 0.005517,
                 0.005517, 0.005517, 0.005517, 0.005517, 0.014996, 0.006718)
  phi6.list[match(T, Period.list)]
}

phi7.cy <- function(T) {
  Period.list <- periods.cy()
  phi7.list <- c(580, 580, 580, 579.9, 579.9, 579.6, 579.2, 577.2, 573.9,
                 568.5, 560.5, 540, 512.9, 441.9, 391.8, 348.1, 332.5, 324.1,
                 321.7, 320.9, 320.3, 320.1, 580, 459)
  phi7.list[match(T, Period.list)]
}

phi8.cy <- function(T) {
  Period.list <- periods.cy()
  phi8.list <- c(0.07, 0.0699, 0.0701, 0.0702, 0.0701, 0.0686, 0.0646, 0.0494,
                 -0.0019, -0.0479, -0.0756, -0.096, -0.0998, -0.0765, -0.0412, 0.014,
                 0.0544, 0.1232, 0.1859, 0.2295, 0.266, 0.2682, 0.07, 0.1138)
  phi8.list[match(T, Period.list)]
}


# 1d. Coefficients for standard deviation term
#     Table 4 in Chiou and Youngs (2008)

tau1.cy <- function(T) {
  Period.list <- periods.cy()
  tau1.list <- c(0.3437, 0.3471, 0.3603, 0.3718, 0.3848, 0.3878, 0.3835, 0.3719,
                 0.3601, 0.3522, 0.3438, 0.3351, 0.3353, 0.3429, 0.3577, 0.3769,
                 0.4023, 0.4406, 0.4784, 0.5074, 0.5328, 0.5542, 0.3437, 0.2539)
  tau1.list[match(T, Period.list)]
}

tau2.cy <- function(T) {
  Period.list <- periods.cy()
  tau2.list <- c(0.2637, 0.2671, 0.2803, 0.2918, 0.3048, 0.3129, 0.3152, 0.3128,
                 0.3076, 0.3047, 0.3005, 0.2984, 0.3036, 0.3205, 0.3419, 0.3703,
                 0.4023, 0.4406, 0.4784, 0.5074, 0.5328, 0.5542, 0.2637, 0.2381)
  tau2.list[match(T, Period.list)]
}

sigma1.cy <- function(T) {
  Period.list <- periods.cy()
  sigma1.list <- c(0.4458, 0.4458, 0.4535, 0.4589, 0.4630, 0.4702, 0.4747, 0.4798,
                   0.4816, 0.4815, 0.4801, 0.4758, 0.4710, 0.4621, 0.4581, 0.4493,
                   0.4459, 0.4433, 0.4424, 0.4420, 0.4416, 0.4414, 0.4458, 0.4496)
  sigma1.list[match(T, Period.list)]
}

sigma2.cy <- function(T) {
  Period.list <- periods.cy()
  sigma2.list <- c(0.3459, 0.3459, 0.3537, 0.3592, 0.3635, 0.3713, 0.3769, 0.3847,
                   0.3902, 0.3946, 0.3981, 0.4036, 0.4079, 0.4157, 0.4213, 0.4213,
                   0.4213, 0.4213, 0.4213, 0.4213, 0.4213, 0.4213, 0.3459, 0.3554)
  sigma2.list[match(T, Period.list)]
}

sigma3.cy <- function(T) {
  Period.list <- periods.cy()
  sigma3.list <- c(0.8000, 0.8000, 0.8000, 0.8000, 0.8000, 0.8000, 0.8000, 0.8000,
                   0.8000, 0.7999, 0.7997, 0.7988, 0.7966, 0.7792, 0.7504, 0.7136,
                   0.7035, 0.7006, 0.7001, 0.7000, 0.7000, 0.7000, 0.8000, 0.7504)
  sigma3.list[match(T, Period.list)]
}

sigma4.cy <- function(T) {
  Period.list <- periods.cy()
  sigma4.list <- c(0.0663, 0.0663, 0.0663, 0.0663, 0.0663, 0.0663, 0.0663, 0.0612,
                   0.0530, 0.0457, 0.0398, 0.0312, 0.0255, 0.0175, 0.0133, 0.0090,
                   0.0068, 0.0045, 0.0034, 0.0027, 0.0018, 0.0014, 0.0663, 0.0133)
  sigma4.list[match(T, Period.list)]
}




# 2. NECESSARY FUNCTION FOR CALCULATING MEDIAN GROUND MOTION
#    Ground motion on reference rock site
LnYref.cy <-  function(M, Rrup, Rjb, Rx, Ztor, dip, Frv, Fnm, Fhw, AS, T) {

  # Load period-independent constants
  coefs.list <- coefs.cy()
  for(i in 1:length(coefs.list))
    eval(parse(text = coefs.list[i]))
  
  # Style-of-faulting term (Line 1 of Eqn 13a)
  F.flt <- c1.cy(T) + (c1a.cy(T)*Frv + c1b.cy(T)*Fnm + c7.cy(T)*(Ztor - 4))*(1 - AS) +
    (c10.cy(T) + c7a.cy(T)*(Ztor - 4))*AS

  # Magnitude term (Line 2 of Eqn 13a)
  F.mag <- c2*(M - 6) +
    ((c2 - c3)/cn.cy(T))*log(1 + exp(cn.cy(T)*(cM.cy(T) - M)))

  # Distance term (Lines 3-5 of Eqn 13a)
  F.dist <- c4*log(Rrup + c5.cy(T)*cosh(c6.cy(T)*max((M - cHM), 0))) +
     (c4a - c4)*log(sqrt(Rrup^2 + cRB^2)) +
       (cGamma1.cy(T) + cGamma2.cy(T) / (cosh(max((M - cGamma3), 0))))*Rrup

  # Hanging wall term (Line 6 of Eqn 13a)
  F.hng <- c9.cy(T) * Fhw * tanh(Rx*(cos(dip*pi/180)^2)/c9a.cy(T)) *
    (1 - sqrt(Rjb^2 + Ztor^2) / (Rrup + 0.001))

  return(F.flt + F.mag + F.dist + F.hng)
}



# 3. NECESSARY FUNCTIONS FOR CALCULATING STANDARD DEVIATION TERM

# 3a. Nonlinear site response term (Eqn 21)
NL.cy <- function(Yref, Vs30, T) {
 
  # Coefficient b (Eqn 10)
  b <- phi2.cy(T) * (exp(phi3.cy(T)*(min(Vs30,1130) - 360)) - exp(phi3.cy(T)*(1130-360)))
  
  # Coefficient c (Eqn 10)
  c <- phi4.cy(T)
  
  # Nonlinear term (Eqn 21)
  return(b*Yref / (Yref + c))
}


# 3b. Sigma, intra-event standard deviation
Sigma.cy <- function(M, Yref, Vs30, VsFlag, AS, T){

  # Nonlinear site response term (Eqn 21)
  NL <- NL.cy(Yref, Vs30, T)

  # Flags for measurement of Vs30
  if(VsFlag == 0){   # Inferred Vs30
    Finferred <- 1
    Fmeasured <- 0
  } else{
    if(VsFlag == 1){   # Measured Vs30
      Finferred <- 0
      Fmeasured <- 1
    }
  }

  # Calculation of sigma (Eqn 20)
  return((sigma1.cy(T) + 0.5*(sigma2.cy(T) - sigma1.cy(T))*(min(max(M,5),7)-5) +
          sigma4.cy(T)*AS) * sqrt(sigma3.cy(T)*Finferred + 0.7*Fmeasured + (1 + NL)^2))
}


# 3c. Tau, inter-event standard deviation (Eqn 19)
Tau.cy <- function(M, T) {
  tau1.cy(T) + 0.5*(tau2.cy(T) - tau1.cy(T)) * (min(max(M,5),7)-5)
}


# 3d. SigmaTot, total standard deviation (Eqn 21)
SigmaTot.cy <- function(M, Yref, Vs30, VsFlag, AS, T) {

  # Intra-event standard deviation
  Sigma <- Sigma.cy(M, Yref, Vs30, VsFlag, AS, T)

  # Inter-event standard deviation
  Tau <- Tau.cy(M, T)

  # Nonlinear site response term
  NL <- NL.cy(Yref, Vs30, T)

  # Total standard deviation
  return(sqrt(((1 + NL)^2)*(Tau^2) + Sigma^2))
}




# 4. MEDIAN GROUND MOTION CALCULATION OF Sa, PGA, and PGV (Eqn 13b)
#    Site Response Function (Eqn 13b)
SaMedian.cy <- function(M, Rjb, Rrup, Rx, Ztor, Frv, Fnm, AS, Fhw, dip, Vs30, Z1.0, T) {
 
  # Ground motion on rock site
  LnYref <- LnYref.cy(M, Rrup, Rjb, Rx, Ztor, dip, Frv, Fnm, Fhw, AS, T)
  Yref <- exp(LnYref)

  # Site response term (Lines 1-2 of Eqn 13b)
  F.site <- phi1.cy(T) * min(log(Vs30/1130), 0) +
    phi2.cy(T) * (exp(phi3.cy(T)*(min(Vs30,1130) - 360)) - exp(phi3.cy(T)*(1130-360))) *
      log((Yref + phi4.cy(T)) / phi4.cy(T))

  # Sediment depth term (Line 3 of Eqn 13b)
  # Slightly modified to prevent COSH error per Chiou
  F.sed <- phi5.cy(T) * (1 - 1 /(cosh(phi6.cy(T) * max(0, (Z1.0 - phi7.cy(T)))))) +
    phi8.cy(T) / cosh(0.15*min(max(0, (Z1.0 - 15)), 300))
  
  return(exp(LnYref + F.site + F.sed))
}




# 5. FINAL FUNCTION FOR CY08 GROUND MOTION CALCULATIONS

Sa.cy <- function(M, Rjb, Vs30, VsFlag, epsilon, T, Rrup = NA, Rx = NA,
                  dip = NA, W = NA, Ztor = NA, Z1.0 = NA, rake = NA, Frv = NA,
                  Fnm = NA, Fhw = NA, azimuth = NA, Zhyp = NA, AS = 0){

  
  # If T is a vector, perform calculation for each of the elements
  if(length(T) > 1 ) {
    return(sapply(T, Sa.cy, M = M, Rjb = Rjb, Vs30 = Vs30, VsFlag = VsFlag,
                  epsilon = epsilon, Rrup = Rrup, Rx = Rx, dip = dip, W = W,
                  Ztor = Ztor, Z1.0 = Z1.0, rake = rake, Frv = Frv, Fnm = Fnm,
                  Fhw = Fhw, azimuth = azimuth, Zhyp = Zhyp, AS = AS))

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
    if(is.na(AS) == TRUE | !(AS == 1 | AS == 0))
      stop("AS must be either 1 (for aftershocks) or 0 (for mainshocks, the default)")
    
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
        if(azimuth >= 0 & azimuth <= 180)
          Fhw <- 1
        else
          Fhw <- 0
      # Calculate from Rx if azimuth is not provided
      } else{
        if(is.na(Rx) == FALSE){
          if(Rx >= 0)
            Fhw <- 1
          else
            Fhw <- 0
        }
      }

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
      Z1.0 <- Z1.calc.cy(Vs30)

  
    # 5C. CALCULATE GROUND MOTION PARAMETER
  
    # Is interpolation necessary?
    interp <- getPeriod(T, "CY08")$interp

    # If interpolation is not necessary, compute Sa
    if(interp == FALSE){
      LnSaMedian <- log(SaMedian.cy(M, Rjb, Rrup, Rx, Ztor, Frv, Fnm,
                                    AS, Fhw, dip, Vs30, Z1.0, T))
      Yref <- exp(LnYref.cy(M, Rrup, Rjb, Rx, Ztor, dip, Frv, Fnm, Fhw, AS, T))
      epsilon.sigmaTot <- epsilon * SigmaTot.cy(M, Yref, Vs30, VsFlag, AS, T)
      LnSa <- LnSaMedian + epsilon.sigmaTot
      return(exp(LnSa))
    } else{

    # If interpolation is necessary, compute Sa
      if(interp == TRUE){
        T1 <- getPeriod(T, "CY08")$lower
        T2 <- getPeriod(T, "CY08")$upper

        # Calculation for T1
        LnSaMedian.T1 <- log(SaMedian.cy(M, Rjb, Rrup, Rx, Ztor, Frv, Fnm,
                                         AS, Fhw, dip, Vs30, Z1.0, T1))
        Yref.T1 <- exp(LnYref.cy(M, Rrup, Rjb, Rx, Ztor, dip, Frv, Fnm, Fhw, AS, T1))
        epsilon.sigmaTot.T1 <- epsilon * SigmaTot.cy(M, Yref.T1, Vs30, VsFlag, AS, T1)
        LnSaT1 <- LnSaMedian.T1 + epsilon.sigmaTot.T1
      
        # Calculation for T2
        LnSaMedian.T2 <- log(SaMedian.cy(M, Rjb, Rrup, Rx, Ztor, Frv, Fnm,
                                         AS, Fhw, dip, Vs30, Z1.0, T2))
        Yref.T2 <- exp(LnYref.cy(M, Rrup, Rjb, Rx, Ztor, dip, Frv, Fnm, Fhw, AS, T2))
        epsilon.sigmaTot.T2 <- epsilon * SigmaTot.cy(M, Yref.T2, Vs30, VsFlag, AS, T2)
        LnSaT2 <- LnSaMedian.T2 + epsilon.sigmaTot.T2
      
        # Interpolated value
        LnSa <- interpolate(log(T), log(T1), log(T2), LnSaT1, LnSaT2)
        return(exp(LnSa))
      }
    }
  }
}

