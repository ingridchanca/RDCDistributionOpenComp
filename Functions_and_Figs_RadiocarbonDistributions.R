# Reproducinng the figures (on main text and supp mat) of the manuscript
# on "Probability distributions of radiocarbon in open linear 
# compartmental systems at steady-state" through the algorithm
# presented in the main text.
# radiocarbon = 14C = RDC

## Packages used in this script
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("SoilR")) install.packages("SoilR")
if (!require("FME")) install.packages("FME")
if (!require("Hmisc")) install.packages("Hmisc")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("grid")) install.packages("grid")
if (!require("gridExtra")) install.packages("gridExtra")
if (!require("cowplot")) install.packages("cowplot")
if (!require("ggpubr")) install.packages("ggpubr")

## Libraries
library(tidyverse)
library(SoilR) #for calculations of distributions of ages, transit time and 14C
library(FME) #for the calculation of distrib (works with SoilR)
library(Hmisc) #statistics (use for variance of theoretical 14C)
library(ggplot2) #for plots
library(grid) #arranging plots
library(gridExtra) #arranging plots
library(cowplot) #arranging plots
library(ggpubr) #for plots, e.g. annotating figs

## directory - these are relative paths depending on the location of this repository
dirname(rstudioapi::getSourceEditorContext()$path)

## plot theme used in the manuscript
theme_set(theme_classic(base_size = 12))

## Functions
### PoolRDC - mass density distribution of 14C in the pools
PoolRDC <- function(A,
                    u,
                    nyears,
                    yrobs,
                    C14atm,
                    C14YearAD,
                    h = 0.1) {
  year.ad <- seq(yrobs, (yrobs - nyears), by = -h)
  age <- seq((yrobs - yrobs), nyears, by = h)
  xss <-  -1 * solve(A) %*% u
  PA <- systemAge(A, u, a = age)$poolAgeDensity
  PAMD <- sapply(1:ncol(A), function(x)
    PA[, x] * xss[x])
  splineC14curve <- splinefun(x = C14YearAD, y = C14atm)
  newResC14curve <- splineC14curve(year.ad)
  D14C.decay <-
    ((((
      newResC14curve / 1000
    ) + 1) * exp((-1 / 8267) * age)) - 1) * 1000
  return(data.frame(
    YearAD = year.ad,
    D14C = D14C.decay,
    Mass = PAMD * h
  ))
}
### SystemRDC - mass density distribution of 14C in the whole system
SystemRDC <- function(A,
                      u,
                      nyears,
                      yrobs,
                      C14atm,
                      C14YearAD,
                      h = 0.1) {
  year.ad <- seq(yrobs, (yrobs - nyears), by = -h)
  age <- seq((yrobs - yrobs), nyears, by = h)
  xss <- -1 * solve(A) %*% u
  SAMD <-
    sum(xss) * (systemAge(A, u, a = age)$systemAgeDensity) #system age mass density
  splineC14curve <- splinefun(x = C14YearAD, y = C14atm)
  newResC14curve <- splineC14curve(year.ad)
  D14C.decay <-
    ((((
      newResC14curve / 1000
    ) + 1) * exp((-1 / 8267) * age)) - 1) * 1000
  return(data.frame(
    YearAD = year.ad,
    D14C = D14C.decay,
    Mass = SAMD * h
  ))
}
### TTRDC - mass density distribution of 14C in the outflux; uses transit time dist
TTRDC <- function(A,
                  u,
                  nyears,
                  yrobs,
                  C14atm,
                  C14YearAD,
                  h = 0.1) {
  year.ad <- seq(yrobs, (yrobs - nyears), by = -h)
  tt <- seq((yrobs - yrobs), nyears, by = h)
  xss <- -1 * solve(A) %*% u
  TTAMD <- sum(xss) * (transitTime(A, u, a = tt)$transitTimeDensity)
  splineC14curve <- splinefun(x = C14YearAD, y = C14atm)
  newResC14curve <- splineC14curve(year.ad)
  D14C.decay <-
    ((((
      newResC14curve / 1000
    ) + 1) * exp((-1 / 8267) * tt)) - 1) * 1000
  return(data.frame(
    YearAD = year.ad,
    D14C = D14C.decay,
    Mass = TTAMD * h
  ))
}
### C14hist - it aggregates the mass densities from previous functions based on their Delta14C
C14hist <- function(D14C,
                    Mass,
                    bin) {
  interval <- seq(floor(min(D14C)),
                  ceiling(max(D14C)),
                  by = bin)
  massdf <- data.frame(
    M = Mass,
    C14class = cut(D14C, breaks = interval),
    labels = FALSE
  )
  aggMass <- aggregate(M ~ C14class,
                       massdf,
                       FUN = sum,
                       drop = FALSE)
  mid <- (diff(interval) / 2) + interval[-length(interval)]
  mass_mid <- cbind(aggMass, mid)
  return(data.frame(mass_mid))
}

## Data

### Soil respiration data - Harvard Forest - Gaudinski et al. (2000), Sierra et al. (2012)
HSMdf <- read.csv(file = "./data/HF14C.csv")

### Atmospheric Delta14CO2
#### Tropics -> IntCal20 + Graven et al. (2017), GMD = records for the tropics
#### + forecast RCP8.5 Graven (2015), PNAS
C14_Trpcs <- read.table(file = "./data/C14Graven_Tropics_forecast_IntCal20") # use for Porce model
#### Northern Hemisphere -> IntCal20 + Graven et al. (2017), GMD = records for NH
#### + forecast RCP8.5 Graven (2015), PNAS
C14_NH <- read.table(file = "./data/C14Graven_PNAS_GMD_Int20_2100_-53050") # use for HFS and Emanuel models

### Compartmental models - inputs and transfer coefficients
#### Harvard Forest Soil (HFS) model - Sierra et al. (2012)
poolnamesHFS <- c(
  "Dead Roots",
  "Oi",
  "Oe/a L",
  "Oe/a H",
  "A, LF (> 80 µm)",
  "A, LF (< 80 µm)",
  "Mineral Associated"
) # names of the pools in order x1, x2, x3, x4, x5, x6 and x7

LI = 150 #Litter inputs
RI = 255 #Root inputs
uH = matrix(nrow = 7,
            ncol = 1, 
            c(RI, LI, 0, 0, 0, 0, 0)) # input vector

# compartmental matrix BH
ks = c(kr = 1/6, koi = 1/1.5, koeal = 1/4, koeah = 1/35,
       kA1 = 1/3, kA2 = 1/75, kM = 1/110)
BH = -abs(diag(ks)) # decomposition rates
BH[3, 2] = ks[2] * (98/(2.7 + 98 + 51)) # transfer from pool 2 to pool 3
BH[4, 3] = ks[3] * (4/(94 + 4)) # transfer from pool 3 to pool 4
BH[6, 5] = ks[5] * (24/(6 + 24)) # transfer from pool 5 to pool 6
BH[7, 6] = ks[6] * (2.4/(22 + 2.4)) # transfer from pool 6 to pool 7
BH[7, 2] = ks[2] * (2.7/(2.7 + 98 + 51)) # transfer from pool 2 to pool 7
BH[4, 1] = ks[1] * (35/(35 + 190 + 30)) # transfer from pool 1 to pool 4
BH[5, 1] = ks[1] * (30/(35 + 190 + 30)) # transfer from pool 1 to pool 5

#### Emanuel model - Emanuel et al. (1981)
poolnamesEm <- c(
  "Non-woody Tree Parts",
  "Woody Tree Parts",
  "Ground Vegetation",
  "Detritus/Decomposers",
  "Active Soil Carbon"
) # names of the pools in order x1, x2, x3, x4 and x5

uE <- matrix(c(77, 0, 36, 0, 0), 5, 1) # input vector

BE <- matrix(
  c(
    -77 / 37, 0, 0, 0, 0,
    31 / 37, -31 / 452, 0, 0, 0,
    0, 0, -36 / 69, 0, 0,
    21 / 37, 15 / 452, 12 / 69, -48 / 81, 0,
    0, 2 / 452, 6 / 69, 3 / 81, -11 / 1121
  ),
  nrow = 5,
  ncol = 5,
  byrow = TRUE
) # compartmental matrix - off-diagonal: transfers between pools;
# diagonal: decomposition rates

#### Porce model - Sierra et al. (2021)
poolnamesP7 <- c(
  "Foliage",
  "Wood",
  "Fine roots",
  "Coarse roots",
  "Fine litter",
  "Coarse woody debris",
  "Soil carbon (0 - 30 cm)"
) # names of the pools in order x1, x2, x3, x4, x5, x6 and x7

GPP <- 23.74244 # GPP input into foliage

uP7 <- matrix(c(GPP, rep(0, 6)), 7, 1) # input vector

# compartmental matrix
BP7 <- diag(
  -1 * c(
    2.97764528,
    0.034705447,
    0.02729699,
    2.221192e-02,
    2.5939474,
    0.5190435,
    0.02350812
  )
) # decomposition rates

BP7[2, 1] = 0.47086499 # transfer from pool 1 to pool 2
BP7[3, 1] = 0.02761772 # transfer from pool 1 to pool 3
BP7[4, 1] = 0.09230466 # transfer from pool 1 to pool 4
BP7[5, 1] = 0.74850986 # transfer from pool 1 to pool 5
BP7[5, 3] = 0.02720847 # transfer from pool 3 to pool 5
BP7[6, 2] = 0.008626976 # transfer from pool 2 to pool 6
BP7[6, 4] = 2.034443e-05 # transfer from pool 4 to pool 6
BP7[7, 5] = 0.6634699 # transfer from pool 5 to pool 7
BP7[7, 6] = 0.5127047 # transfer from pool 6 to pool 7

## Calculations
age <- seq(0, 1000, by = 0.1) # time span for age and transit time distributions
### HFS model
#### Age and Transit Time distributions
SAH <- systemAge(A = BH, u = uH, a = age)
SAHpool <- data.frame(SAH$poolAgeDensity) # age distributions of individual pools
SAHsys <- data.frame(SAH$systemAgeDensity) # age distribution of the whole system
TTH <- transitTime(A = BH, u = uH, a = age)
TTHdens <- data.frame(TTH$transitTimeDensity) # transit time distribution - age of outflux

#### Radiocarbon distributions
##### Year of observation = AD 1965
HPoolRDC1965 <- PoolRDC(
  A = BH,
  u = uH,
  nyears = 1000,
  yrobs = 1965,
  C14atm = C14_NH$Delta.14C,
  C14YearAD = C14_NH$Year.AD,
  h = 0.1
) # for pools

HP1hist1965 <- C14hist(D14C = HPoolRDC1965$D14C,
                       Mass = HPoolRDC1965$Mass.1,
                       bin = 40) # aggregation for pool 1
HP2hist1965 <- C14hist(D14C = HPoolRDC1965$D14C,
                       Mass = HPoolRDC1965$Mass.2,
                       bin = 40) # aggregation for pool 2
HP3hist1965 <- C14hist(D14C = HPoolRDC1965$D14C,
                       Mass = HPoolRDC1965$Mass.3,
                       bin = 40) # aggregation for pool 3
HP4hist1965 <- C14hist(D14C = HPoolRDC1965$D14C,
                       Mass = HPoolRDC1965$Mass.4,
                       bin = 40) # aggregation for pool 4
HP5hist1965 <- C14hist(D14C = HPoolRDC1965$D14C,
                       Mass = HPoolRDC1965$Mass.5,
                       bin = 40) # aggregation for pool 5
HP6hist1965 <- C14hist(D14C = HPoolRDC1965$D14C,
                       Mass = HPoolRDC1965$Mass.6,
                       bin = 40) # aggregation for pool 6
HP7hist1965 <- C14hist(D14C = HPoolRDC1965$D14C,
                       Mass = HPoolRDC1965$Mass.7,
                       bin = 40) # aggregation for pool 7

HP1.1965 <-
  weighted.mean(x = HPoolRDC1965$D14C,
                w = HPoolRDC1965$Mass.1) # expected value for pool 1
HP2.1965 <-
  weighted.mean(x = HPoolRDC1965$D14C,
                w = HPoolRDC1965$Mass.2) # expected value for pool 2
HP3.1965 <-
  weighted.mean(x = HPoolRDC1965$D14C,
                w = HPoolRDC1965$Mass.3) # expected value for pool 3
HP4.1965 <-
  weighted.mean(x = HPoolRDC1965$D14C,
                w = HPoolRDC1965$Mass.4) # expected value for pool 4
HP5.1965 <-
  weighted.mean(x = HPoolRDC1965$D14C,
                w = HPoolRDC1965$Mass.5) # expected value for pool 5
HP6.1965 <-
  weighted.mean(x = HPoolRDC1965$D14C,
                w = HPoolRDC1965$Mass.6) # expected value for pool 6
HP7.1965 <-
  weighted.mean(x = HPoolRDC1965$D14C,
                w = HPoolRDC1965$Mass.7) # expected value for pool 7

RDCexpPoolsH.1965 <- data.frame(HP1.1965,
                                HP2.1965,
                                HP3.1965,
                                HP4.1965,
                                HP5.1965,
                                HP6.1965,
                                HP7.1965,
                                row.names = "Δ14C")
colnames(RDCexpPoolsH.1965) <- c(
  poolnamesHFS[1],
  poolnamesHFS[2],
  poolnamesHFS[3],
  poolnamesHFS[4],
  poolnamesHFS[5],
  poolnamesHFS[6],
  poolnamesHFS[7]
)
RDCexpPoolsH.1965 # summary of expected values for pools in a table

HSysRDC1965 <- SystemRDC(
  A = BH,
  u = uH,
  nyears = 1000,
  yrobs = 1965,
  C14atm = C14_NH$Delta.14C,
  C14YearAD = C14_NH$Year.AD,
  h = 0.1
) # for the whole system

HSyshist1965 <- C14hist(D14C = HSysRDC1965$D14C,
                        Mass = HSysRDC1965$Mass,
                        bin = 40) # aggregation for the whole system

HSys.1965 <-
  weighted.mean(x = HSysRDC1965$D14C,
                w = HSysRDC1965$Mass) # expected value for the whole system
var.HSys1965 <-
  wtd.var(x = HSysRDC1965$D14C, 
          weights = HSysRDC1965$Mass)
sd.HSys1965 <- sqrt(var.HSys1965) # standard deviation of the distribution

RDCexpSysH.1965 <- data.frame(HSys.1965,
                             sd.HSys1965,
                             row.names = "System")
colnames(RDCexpSysH.1965) <- c("Δ14C",
                              "sd")
RDCexpSysH.1965 # summary of expected value and sd of system in a table

HTTRDC1965 <- TTRDC(
  A = BH,
  u = uH,
  nyears = 1000,
  yrobs = 1965,
  C14atm = C14_NH$Delta.14C,
  C14YearAD = C14_NH$Year.AD,
  h = 0.1
) # for the outflux, i.e. based on the transit time distribution

HTThist1965 <- C14hist(D14C = HTTRDC1965$D14C,
                       Mass = HTTRDC1965$Mass,
                       bin = 40) # aggregation for the outflux

HTT.1965 <-
  weighted.mean(x = HTTRDC1965$D14C,
                w = HTTRDC1965$Mass) # expected value of the outflux
var.HTT1965 <-
  wtd.var(x = HTTRDC1965$D14C, 
          weights = HTTRDC1965$Mass)
sd.HTT1965 <- sqrt(var.HTT1965) # standard deviation of the distribution

RDCexpTTH.1965 <- data.frame(HTT.1965,
                             sd.HTT1965,
                             row.names = "Outflux")
colnames(RDCexpTTH.1965) <- c("Δ14C",
                              "sd")
RDCexpTTH.1965 # summary of expected value and sd of outflux in a table

##### Year of observation = AD 1996
HTTRDC1996 <- TTRDC(
  A = BH,
  u = uH,
  nyears = 1000,
  yrobs = 1996,
  C14atm = C14_NH$Delta.14C,
  C14YearAD = C14_NH$Year.AD,
  h = 0.1
) # for outflux, i.e. based on transit time distribution

HTThist1996 <- C14hist(D14C = HTTRDC1996$D14C,
                       Mass = HTTRDC1996$Mass,
                       bin = 10) # aggregation for the outflux

HTT.1996 <-
  weighted.mean(x = HTTRDC1996$D14C,
                w = HTTRDC1996$Mass) # expected value of outflux - year 1996
var.HTT1996 <-
  wtd.var(x = HTTRDC1996$D14C, 
          weights = HTTRDC1996$Mass)
sd.HTT1996 <- sqrt(var.HTT1996) # standard deviation of the distribution

##### Year of observation = AD 1998
HTTRDC1998 <- TTRDC(
  A = BH,
  u = uH,
  nyears = 1000,
  yrobs = 1998,
  C14atm = C14_NH$Delta.14C,
  C14YearAD = C14_NH$Year.AD,
  h = 0.1
) # for outflux, i.e. based on transit time distribution

HTThist1998 <- C14hist(D14C = HTTRDC1998$D14C,
                       Mass = HTTRDC1998$Mass,
                       bin = 10) # aggregation for the outflux

HTT.1998 <-
  weighted.mean(x = HTTRDC1998$D14C,
                w = HTTRDC1998$Mass) # expected value of outflux - year 1998
var.HTT1998 <-
  wtd.var(x = HTTRDC1998$D14C, 
          weights = HTTRDC1998$Mass)
sd.HTT1998 <- sqrt(var.HTT1998) # standard deviation of the distribution

##### Year of observation = AD 2002
HTTRDC2002 <- TTRDC(
  A = BH,
  u = uH,
  nyears = 1000,
  yrobs = 2002,
  C14atm = C14_NH$Delta.14C,
  C14YearAD = C14_NH$Year.AD,
  h = 0.1
) # for outflux, i.e. based on transit time distribution

HTThist2002 <- C14hist(D14C = HTTRDC2002$D14C,
                       Mass = HTTRDC2002$Mass,
                       bin = 10) # aggregation for the outflux
HTT.2002 <-
  weighted.mean(x = HTTRDC2002$D14C,
                w = HTTRDC2002$Mass) # expected value of outflux - year 2002
var.HTT2002 <-
  wtd.var(x = HTTRDC2002$D14C, 
          weights = HTTRDC2002$Mass)
sd.HTT2002 <- sqrt(var.HTT2002) # standard deviation of the distribution

##### Year of observation = AD 2008
HTTRDC2008 <- TTRDC(
  A = BH,
  u = uH,
  nyears = 1000,
  yrobs = 2008,
  C14atm = C14_NH$Delta.14C,
  C14YearAD = C14_NH$Year.AD,
  h = 0.1
) # for outflux, i.e. based on transit time distribution

HTThist2008 <- C14hist(D14C = HTTRDC2008$D14C,
                       Mass = HTTRDC2008$Mass,
                       bin = 10) # aggregation for the outflux

HTT.2008 <-
  weighted.mean(x = HTTRDC2008$D14C,
                w = HTTRDC2008$Mass) # expected value of outflux - year 2008
var.HTT2008 <-
  wtd.var(x = HTTRDC2008$D14C, 
          weights = HTTRDC2008$Mass)
sd.HTT2008 <- sqrt(var.HTT2008) # standard deviation of the distribution

##### Year of observation = AD 2027
HPoolRDC2027 <- PoolRDC(
  A = BH,
  u = uH,
  nyears = 1000,
  yrobs = 2027,
  C14atm = C14_NH$Delta.14C,
  C14YearAD = C14_NH$Year.AD,
  h = 0.1
) # for pools

HP1hist2027 <- C14hist(D14C = HPoolRDC2027$D14C,
                       Mass = HPoolRDC2027$Mass.1,
                       bin = 10) # aggregation for pool 1
HP2hist2027 <- C14hist(D14C = HPoolRDC2027$D14C,
                       Mass = HPoolRDC2027$Mass.2,
                       bin = 10) # aggregation for pool 2
HP3hist2027 <- C14hist(D14C = HPoolRDC2027$D14C,
                       Mass = HPoolRDC2027$Mass.3,
                       bin = 10) # aggregation for pool 3
HP4hist2027 <- C14hist(D14C = HPoolRDC2027$D14C,
                       Mass = HPoolRDC2027$Mass.4,
                       bin = 10) # aggregation for pool 4
HP5hist2027 <- C14hist(D14C = HPoolRDC2027$D14C,
                       Mass = HPoolRDC2027$Mass.5,
                       bin = 10) # aggregation for pool 5
HP6hist2027 <- C14hist(D14C = HPoolRDC2027$D14C,
                       Mass = HPoolRDC2027$Mass.6,
                       bin = 10) # aggregation for pool 6
HP7hist2027 <- C14hist(D14C = HPoolRDC2027$D14C,
                       Mass = HPoolRDC2027$Mass.7,
                       bin = 10) # aggregation for pool 7

HP1.2027 <-
  weighted.mean(x = HPoolRDC2027$D14C,
                w = HPoolRDC2027$Mass.1) # expected value of pool 1
HP2.2027 <-
  weighted.mean(x = HPoolRDC2027$D14C,
                w = HPoolRDC2027$Mass.2) # expected value of pool 2
HP3.2027 <-
  weighted.mean(x = HPoolRDC2027$D14C,
                w = HPoolRDC2027$Mass.3) # expected value of pool 3
HP4.2027 <-
  weighted.mean(x = HPoolRDC2027$D14C,
                w = HPoolRDC2027$Mass.4) # expected value of pool 4
HP5.2027 <-
  weighted.mean(x = HPoolRDC2027$D14C,
                w = HPoolRDC2027$Mass.5) # expected value of pool 5
HP6.2027 <-
  weighted.mean(x = HPoolRDC2027$D14C,
                w = HPoolRDC2027$Mass.6) # expected value of pool 6
HP7.2027 <-
  weighted.mean(x = HPoolRDC2027$D14C,
                w = HPoolRDC2027$Mass.7) # expected value of pool 7

RDCexpPoolsH.2027 <- data.frame(HP1.2027,
                                HP2.2027,
                                HP3.2027,
                                HP4.2027,
                                HP5.2027,
                                HP6.2027,
                                HP7.2027,
                                row.names = "Δ14C")
colnames(RDCexpPoolsH.2027) <- c(
  poolnamesHFS[1],
  poolnamesHFS[2],
  poolnamesHFS[3],
  poolnamesHFS[4],
  poolnamesHFS[5],
  poolnamesHFS[6],
  poolnamesHFS[7]
)

RDCexpPoolsH.2027 # summary of expected values of pools in a table

HSysRDC2027 <- SystemRDC(
  A = BH,
  u = uH,
  nyears = 1000,
  yrobs = 2027,
  C14atm = C14_NH$Delta.14C,
  C14YearAD = C14_NH$Year.AD,
  h = 0.1
) # for the whole system

HSyshist2027 <- C14hist(D14C = HSysRDC2027$D14C,
                        Mass = HSysRDC2027$Mass,
                        bin = 10) # aggregation for system
HSys.2027 <-
  weighted.mean(x = HSysRDC2027$D14C,
                w = HSysRDC2027$Mass) # expected value of system
var.HSys2027 <-
  wtd.var(x = HSysRDC2027$D14C, 
          weights = HSysRDC2027$Mass)
sd.HSys2027 <- sqrt(var.HSys2027) # standard deviation of the distribution

HTTRDC2027 <- TTRDC(
  A = BH,
  u = uH,
  nyears = 1000,
  yrobs = 2027,
  C14atm = C14_NH$Delta.14C,
  C14YearAD = C14_NH$Year.AD,
  h = 0.1
) # for outflux, i.e. based on transit time distribution

HTThist2027 <- C14hist(D14C = HTTRDC2027$D14C,
                       Mass = HTTRDC2027$Mass,
                       bin = 10) # aggregation of outflux

HTThist2027.p <- C14hist(D14C = HTTRDC2027$D14C,
                         Mass = HTTRDC2027$Mass,
                         bin = 40) # with bin size b=40

HTT.2027 <-
  weighted.mean(x = HTTRDC2027$D14C,
                w = HTTRDC2027$Mass) # expected value of outflux
var.HTT2027 <-
  wtd.var(x = HTTRDC2027$D14C, 
          weights = HTTRDC2027$Mass)
sd.HTT2027 <- sqrt(var.HTT2027) # standard deviation of the distribution

##### Year of observation = AD 2100
HPoolRDC2100 <- PoolRDC(
  A = BH,
  u = uH,
  nyears = 1000,
  yrobs = 2100,
  C14atm = C14_NH$Delta.14C,
  C14YearAD = C14_NH$Year.AD,
  h = 0.1
) # for pools

HP1hist2100 <- C14hist(D14C = HPoolRDC2100$D14C,
                       Mass = HPoolRDC2100$Mass.1,
                       bin = 10) # aggregation for pool 1
HP2hist2100 <- C14hist(D14C = HPoolRDC2100$D14C,
                       Mass = HPoolRDC2100$Mass.2,
                       bin = 10) # aggregation for pool 2
HP3hist2100 <- C14hist(D14C = HPoolRDC2100$D14C,
                       Mass = HPoolRDC2100$Mass.3,
                       bin = 10) # aggregation for pool 3
HP4hist2100 <- C14hist(D14C = HPoolRDC2100$D14C,
                       Mass = HPoolRDC2100$Mass.4,
                       bin = 10) # aggregation for pool 4
HP5hist2100 <- C14hist(D14C = HPoolRDC2100$D14C,
                       Mass = HPoolRDC2100$Mass.5,
                       bin = 10) # aggregation for pool 5
HP6hist2100 <- C14hist(D14C = HPoolRDC2100$D14C,
                       Mass = HPoolRDC2100$Mass.6,
                       bin = 10) # aggregation for pool 6
HP7hist2100 <- C14hist(D14C = HPoolRDC2100$D14C,
                       Mass = HPoolRDC2100$Mass.7,
                       bin = 10) # aggregation for pool 7

HP1.2100 <-
  weighted.mean(x = HPoolRDC2100$D14C,
                w = HPoolRDC2100$Mass.1) # expected value of pool 1
HP2.2100 <-
  weighted.mean(x = HPoolRDC2100$D14C,
                w = HPoolRDC2100$Mass.2) # expected value of pool 2
HP3.2100 <-
  weighted.mean(x = HPoolRDC2100$D14C,
                w = HPoolRDC2100$Mass.3) # expected value of pool 3
HP4.2100 <-
  weighted.mean(x = HPoolRDC2100$D14C,
                w = HPoolRDC2100$Mass.4) # expected value of pool 4
HP5.2100 <-
  weighted.mean(x = HPoolRDC2100$D14C,
                w = HPoolRDC2100$Mass.5) # expected value of pool 5
HP6.2100 <-
  weighted.mean(x = HPoolRDC2100$D14C,
                w = HPoolRDC2100$Mass.6) # expected value of pool 6
HP7.2100 <-
  weighted.mean(x = HPoolRDC2100$D14C,
                w = HPoolRDC2100$Mass.7) # expected value of pool 7

RDCexpPoolsH.2100 <- data.frame(HP1.2100,
                                HP2.2100,
                                HP3.2100,
                                HP4.2100,
                                HP5.2100,
                                HP6.2100,
                                HP7.2100,
                                row.names = "Δ14C")
colnames(RDCexpPoolsH.2100) <- c(
  poolnamesHFS[1],
  poolnamesHFS[2],
  poolnamesHFS[3],
  poolnamesHFS[4],
  poolnamesHFS[5],
  poolnamesHFS[6],
  poolnamesHFS[7]
)
RDCexpPoolsH.2100 # summary of expected values of pools in a table

HSysRDC2100 <- SystemRDC(
  A = BH,
  u = uH,
  nyears = 1000,
  yrobs = 2100,
  C14atm = C14_NH$Delta.14C,
  C14YearAD = C14_NH$Year.AD,
  h = 0.1
) # for the whole system 

HSyshist2100 <- C14hist(D14C = HSysRDC2100$D14C,
                        Mass = HSysRDC2100$Mass,
                        bin = 10) # aggregation for system

HSys.2100 <-
  weighted.mean(x = HSysRDC2100$D14C,
                w = HSysRDC2100$Mass) # expected value of system
var.HSys2100 <-
  wtd.var(x = HSysRDC2100$D14C, 
          weights = HSysRDC2100$Mass)
sd.HSys2100 <- sqrt(var.HSys2100) # standard deviation of the distribution

HTTRDC2100 <- TTRDC(
  A = BH,
  u = uH,
  nyears = 1000,
  yrobs = 2100,
  C14atm = C14_NH$Delta.14C,
  C14YearAD = C14_NH$Year.AD,
  h = 0.1
) # for outflux, i.e. based on transit time distribution

HTThist2100 <- C14hist(D14C = HTTRDC2100$D14C,
                       Mass = HTTRDC2100$Mass,
                       bin = 10) # aggregation for outflux

HTThist2100.p <- C14hist(D14C = HTTRDC2100$D14C,
                         Mass = HTTRDC2100$Mass,
                         bin = 40) # with bin size b=40

HTT.2100 <-
  weighted.mean(x = HTTRDC2100$D14C,
                w = HTTRDC2100$Mass) # expected value of outflux
var.HTT2100 <-
  wtd.var(x = HTTRDC2100$D14C, 
          weights = HTTRDC2100$Mass)
sd.HTT2100 <- sqrt(var.HTT2100) # standard deviation of the distribution

### Porce model
#### Age and Transit Time distributions
SAP7 <- systemAge(A = BP7, u = uP7, a = age)
SAP7pool <- data.frame(SAP7$poolAgeDensity) # age distributions of individual pools
SAP7sys <- data.frame(SAP7$systemAgeDensity) # age distribution of the whole system
TTP7 <- transitTime(A = BP7, u = uP7, a = age)
TTP7dens <- data.frame(TTP7$transitTimeDensity) # transit time distribution - age of outflux

#### Radiocarbon distributions
##### Year of observation = AD 1965
P7PoolRDC1965 <- PoolRDC(
  A = BP7,
  u = uP7,
  nyears = 1000,
  yrobs = 1965,
  C14atm = C14_Trpcs$Delta.14C,
  C14YearAD = C14_Trpcs$Year.AD,
  h = 0.1
) # for pools

P7P1hist1965 <- C14hist(D14C = P7PoolRDC1965$D14C,
                        Mass = P7PoolRDC1965$Mass.1,
                        bin = 40) # aggregation for pool 1
P7P2hist1965 <- C14hist(D14C = P7PoolRDC1965$D14C,
                        Mass = P7PoolRDC1965$Mass.2,
                        bin = 40) # aggregation for pool 2
P7P3hist1965 <- C14hist(D14C = P7PoolRDC1965$D14C,
                        Mass = P7PoolRDC1965$Mass.3,
                        bin = 40) # aggregation for pool 3
P7P4hist1965 <- C14hist(D14C = P7PoolRDC1965$D14C,
                        Mass = P7PoolRDC1965$Mass.4,
                        bin = 40) # aggregation for pool 4
P7P5hist1965 <- C14hist(D14C = P7PoolRDC1965$D14C,
                        Mass = P7PoolRDC1965$Mass.5,
                        bin = 40) # aggregation for pool 5
P7P6hist1965 <- C14hist(D14C = P7PoolRDC1965$D14C,
                        Mass = P7PoolRDC1965$Mass.6,
                        bin = 40) # aggregation for pool 6
P7P7hist1965 <- C14hist(D14C = P7PoolRDC1965$D14C,
                        Mass = P7PoolRDC1965$Mass.7,
                        bin = 40) # aggregation for pool 7

P7P1.1965 <-
  weighted.mean(x = P7PoolRDC1965$D14C,
                w = P7PoolRDC1965$Mass.1) # expected value for pool 1
P7P2.1965 <-
  weighted.mean(x = P7PoolRDC1965$D14C,
                w = P7PoolRDC1965$Mass.2) # expected value for pool 2
P7P3.1965 <-
  weighted.mean(x = P7PoolRDC1965$D14C,
                w = P7PoolRDC1965$Mass.3) # expected value for pool 3
P7P4.1965 <-
  weighted.mean(x = P7PoolRDC1965$D14C,
                w = P7PoolRDC1965$Mass.4) # expected value for pool 4
P7P5.1965 <-
  weighted.mean(x = P7PoolRDC1965$D14C,
                w = P7PoolRDC1965$Mass.5) # expected value for pool 5
P7P6.1965 <-
  weighted.mean(x = P7PoolRDC1965$D14C,
                w = P7PoolRDC1965$Mass.6) # expected value for pool 6
P7P7.1965 <-
  weighted.mean(x = P7PoolRDC1965$D14C,
                w = P7PoolRDC1965$Mass.7) # expected value for pool 7

RDCexpPoolsP7.1965 <- data.frame(
  P7P1.1965,
  P7P2.1965,
  P7P3.1965,
  P7P4.1965,
  P7P5.1965,
  P7P6.1965,
  P7P7.1965,
  row.names = "Δ14C"
)
colnames(RDCexpPoolsP7.1965) <- c(
  poolnamesP7[1],
  poolnamesP7[2],
  poolnamesP7[3],
  poolnamesP7[4],
  poolnamesP7[5],
  poolnamesP7[6],
  poolnamesP7[7]
)
RDCexpPoolsP7.1965 # summary of expected values of pools in a table

P7SysRDC1965 <- SystemRDC(
  A = BP7,
  u = uP7,
  nyears = 1000,
  yrobs = 1965,
  C14atm = C14_Trpcs$Delta.14C,
  C14YearAD = C14_Trpcs$Year.AD,
  h = 0.1
) # for the whole system

P7Syshist1965 <- C14hist(D14C = P7SysRDC1965$D14C,
                         Mass = P7SysRDC1965$Mass,
                         bin = 40) # aggregation for system

P7Sys.1965 <-
  weighted.mean(x = P7SysRDC1965$D14C,
                w = P7SysRDC1965$Mass) # expected value of system
var.P7Sys1965 <-
  wtd.var(x = P7SysRDC1965$D14C,
          weights = P7SysRDC1965$Mass)
sd.P7Sys1965 <- sqrt(var.P7Sys1965) # standard deviation of the distribution

RDCexpSysP7.1965 <- data.frame(P7Sys.1965,
                               sd.P7Sys1965,
                               row.names = "Whole System")
colnames(RDCexpSysP7.1965) <- c("Δ14C",
                                "sd")
RDCexpSysP7.1965 # summary of expected value and sd of system in a table

P7TTRDC1965 <- TTRDC(
  A = BP7,
  u = uP7,
  nyears = 1000,
  yrobs = 1965,
  C14atm = C14_Trpcs$Delta.14C,
  C14YearAD = C14_Trpcs$Year.AD,
  h = 0.1
) # for outflux, i.e. based on transit time distribution

P7TThist1965 <- C14hist(D14C = P7TTRDC1965$D14C,
                        Mass = P7TTRDC1965$Mass,
                        bin = 40) # aggregation for outflux

P7TT.1965 <-
  weighted.mean(x = P7TTRDC1965$D14C,
                w = P7TTRDC1965$Mass) # expected value of outflux
var.P7TT1965 <-
  wtd.var(x = P7TTRDC1965$D14C, 
          weights = P7TTRDC1965$Mass)
sd.P7TT1965 <- sqrt(var.P7TT1965) # standard deviation of the distribution

RDCexpTTP7.1965 <- data.frame(P7TT.1965,
                              sd.P7TT1965,
                              row.names = "Outflux")
colnames(RDCexpTTP7.1965) <- c("Δ14C",
                               "sd")
RDCexpTTP7.1965 # summary of expected value and sd of outflux in a table

##### Year of observation = AD 2027
P7PoolRDC2027 <- PoolRDC(
  A = BP7,
  u = uP7,
  nyears = 1000,
  yrobs = 2027,
  C14atm = C14_Trpcs$Delta.14C,
  C14YearAD = C14_Trpcs$Year.AD,
  h = 0.1
) # for pools

P7P1hist2027 <- C14hist(D14C = P7PoolRDC2027$D14C,
                        Mass = P7PoolRDC2027$Mass.1,
                        bin = 10) # aggregation for pool 1
P7P2hist2027 <- C14hist(D14C = P7PoolRDC2027$D14C,
                        Mass = P7PoolRDC2027$Mass.2,
                        bin = 10) # aggregation for pool 2
P7P3hist2027 <- C14hist(D14C = P7PoolRDC2027$D14C,
                        Mass = P7PoolRDC2027$Mass.3,
                        bin = 10) # aggregation for pool 3
P7P4hist2027 <- C14hist(D14C = P7PoolRDC2027$D14C,
                        Mass = P7PoolRDC2027$Mass.4,
                        bin = 10) # aggregation for pool 4
P7P5hist2027 <- C14hist(D14C = P7PoolRDC2027$D14C,
                        Mass = P7PoolRDC2027$Mass.5,
                        bin = 10) # aggregation for pool 5
P7P6hist2027 <- C14hist(D14C = P7PoolRDC2027$D14C,
                        Mass = P7PoolRDC2027$Mass.6,
                        bin = 10) # aggregation for pool 6
P7P7hist2027 <- C14hist(D14C = P7PoolRDC2027$D14C,
                        Mass = P7PoolRDC2027$Mass.7,
                        bin = 10) # aggregation for pool 7

P7P1.2027 <-
  weighted.mean(x = P7PoolRDC2027$D14C,
                w = P7PoolRDC2027$Mass.1) # expected value of pool 1
P7P2.2027 <-
  weighted.mean(x = P7PoolRDC2027$D14C,
                w = P7PoolRDC2027$Mass.2) # expected value of pool 2
P7P3.2027 <-
  weighted.mean(x = P7PoolRDC2027$D14C,
                w = P7PoolRDC2027$Mass.3) # expected value of pool 3
P7P4.2027 <-
  weighted.mean(x = P7PoolRDC2027$D14C,
                w = P7PoolRDC2027$Mass.4) # expected value of pool 4
P7P5.2027 <-
  weighted.mean(x = P7PoolRDC2027$D14C,
                w = P7PoolRDC2027$Mass.5) # expected value of pool 5
P7P6.2027 <-
  weighted.mean(x = P7PoolRDC2027$D14C,
                w = P7PoolRDC2027$Mass.6) # expected value of pool 6
P7P7.2027 <-
  weighted.mean(x = P7PoolRDC2027$D14C,
                w = P7PoolRDC2027$Mass.7) # expected value of pool 7

RDCexpPoolsP7.2027 <- data.frame(
  P7P1.2027,
  P7P2.2027,
  P7P3.2027,
  P7P4.2027,
  P7P5.2027,
  P7P6.2027,
  P7P7.2027,
  row.names = "Δ14C"
)
colnames(RDCexpPoolsP7.2027) <- c(
  poolnamesP7[1],
  poolnamesP7[2],
  poolnamesP7[3],
  poolnamesP7[4],
  poolnamesP7[5],
  poolnamesP7[6],
  poolnamesP7[7]
)
RDCexpPoolsP7.2027 # summary of expected value of pools in a table

P7SysRDC2027 <- SystemRDC(
  A = BP7,
  u = uP7,
  nyears = 1000,
  yrobs = 2027,
  C14atm = C14_Trpcs$Delta.14C,
  C14YearAD = C14_Trpcs$Year.AD,
  h = 0.1
) # for the whole system

P7Syshist2027 <- C14hist(D14C = P7SysRDC2027$D14C,
                         Mass = P7SysRDC2027$Mass,
                         bin = 10) # aggregation of system

P7Sys.2027 <-
  weighted.mean(x = P7SysRDC2027$D14C,
                w = P7SysRDC2027$Mass) # expected value of system
var.P7Sys2027 <-
  wtd.var(x = P7SysRDC2027$D14C,
          weights = P7SysRDC2027$Mass)
sd.P7Sys2027 <- sqrt(var.P7Sys2027) # standard deviation of the distribution

RDCexpSysP7.2027 <- data.frame(P7Sys.2027,
                               sd.P7Sys2027,
                               row.names = "Whole System")
colnames(RDCexpSysP7.2027) <- c("Δ14C",
                                "sd")
RDCexpSysP7.2027 # summary of expected value and sd of system in a table

P7TTRDC2027 <- TTRDC(
  A = BP7,
  u = uP7,
  nyears = 1000,
  yrobs = 2027,
  C14atm = C14_Trpcs$Delta.14C,
  C14YearAD = C14_Trpcs$Year.AD,
  h = 0.1
) # for outflux, i.e. based on transit time distribution

P7TThist2027 <- C14hist(D14C = P7TTRDC2027$D14C,
                        Mass = P7TTRDC2027$Mass,
                        bin = 10) # aggregation for outflux

P7TThist2027.p <- C14hist(D14C = P7TTRDC2027$D14C,
                          Mass = P7TTRDC2027$Mass,
                          bin = 40) # with bin size b=40

P7TT.2027 <-
  weighted.mean(x = P7TTRDC2027$D14C,
                w = P7TTRDC2027$Mass) # expected value of outflux
var.P7TT2027 <-
  wtd.var(x = P7TTRDC2027$D14C, 
          weights = P7TTRDC2027$Mass)
sd.P7TT2027 <- sqrt(var.P7TT2027) # standard deviation of the distribution

RDCexpTTP7.2027 <- data.frame(P7TT.2027,
                              sd.P7TT2027,
                              row.names = "Outflux")
colnames(RDCexpTTP7.2027) <- c("Δ14C",
                               "sd")
RDCexpTTP7.2027 # summary of expected value and sd of outflux in a table

##### Year of observation = AD 2100
P7PoolRDC2100 <- PoolRDC(
  A = BP7,
  u = uP7,
  nyears = 1000,
  yrobs = 2100,
  C14atm = C14_Trpcs$Delta.14C,
  C14YearAD = C14_Trpcs$Year.AD,
  h = 0.1
) # for pools

P7P1hist2100 <- C14hist(D14C = P7PoolRDC2100$D14C,
                        Mass = P7PoolRDC2100$Mass.1,
                        bin = 10) # aggregation for pool 1
P7P2hist2100 <- C14hist(D14C = P7PoolRDC2100$D14C,
                        Mass = P7PoolRDC2100$Mass.2,
                        bin = 10) # aggregation for pool 2
P7P3hist2100 <- C14hist(D14C = P7PoolRDC2100$D14C,
                        Mass = P7PoolRDC2100$Mass.3,
                        bin = 10) # aggregation for pool 3
P7P4hist2100 <- C14hist(D14C = P7PoolRDC2100$D14C,
                        Mass = P7PoolRDC2100$Mass.4,
                        bin = 10) # aggregation for pool 4
P7P5hist2100 <- C14hist(D14C = P7PoolRDC2100$D14C,
                        Mass = P7PoolRDC2100$Mass.5,
                        bin = 10) # aggregation for pool 5
P7P6hist2100 <- C14hist(D14C = P7PoolRDC2100$D14C,
                        Mass = P7PoolRDC2100$Mass.6,
                        bin = 10) # aggregation for pool 6
P7P7hist2100 <- C14hist(D14C = P7PoolRDC2100$D14C,
                        Mass = P7PoolRDC2100$Mass.7,
                        bin = 10) # aggregation for pool 7

P7P1.2100 <-
  weighted.mean(x = P7PoolRDC2100$D14C,
                w = P7PoolRDC2100$Mass.1) # expected value of pool 1
P7P2.2100 <-
  weighted.mean(x = P7PoolRDC2100$D14C,
                w = P7PoolRDC2100$Mass.2) # expected value of pool 2
P7P3.2100 <-
  weighted.mean(x = P7PoolRDC2100$D14C,
                w = P7PoolRDC2100$Mass.3) # expected value of pool 3
P7P4.2100 <-
  weighted.mean(x = P7PoolRDC2100$D14C,
                w = P7PoolRDC2100$Mass.4) # expected value of pool 4
P7P5.2100 <-
  weighted.mean(x = P7PoolRDC2100$D14C,
                w = P7PoolRDC2100$Mass.5) # expected value of pool 5
P7P6.2100 <-
  weighted.mean(x = P7PoolRDC2100$D14C,
                w = P7PoolRDC2100$Mass.6) # expected value of pool 6
P7P7.2100 <-
  weighted.mean(x = P7PoolRDC2100$D14C,
                w = P7PoolRDC2100$Mass.7) # expected value of pool 7

RDCexpPoolsP7.2100 <- data.frame(P7P1.2100, 
                                 P7P2.2100,
                                 P7P3.2100, 
                                 P7P4.2100, 
                                 P7P5.2100,
                                 P7P6.2100,
                                 P7P7.2100,
                                 row.names = "Δ14C")
colnames(RDCexpPoolsP7.2100) <- c(
  poolnamesP7[1],
  poolnamesP7[2],
  poolnamesP7[3],
  poolnamesP7[4],
  poolnamesP7[5],
  poolnamesP7[6],
  poolnamesP7[7]
)
RDCexpPoolsP7.2100 # summary of expected values of pools in a table

P7SysRDC2100 <- SystemRDC(
  A = BP7,
  u = uP7,
  nyears = 1000,
  yrobs = 2100,
  C14atm = C14_Trpcs$Delta.14C,
  C14YearAD = C14_Trpcs$Year.AD,
  h = 0.1
) # for the whole system

P7Syshist2100 <- C14hist(D14C = P7SysRDC2100$D14C,
                         Mass = P7SysRDC2100$Mass,
                         bin = 10) # aggregation for system

P7Sys.2100 <-
  weighted.mean(x = P7SysRDC2100$D14C,
                w = P7SysRDC2100$Mass) # expected value of system
var.P7Sys2100 <-
  wtd.var(x = P7SysRDC2100$D14C,
          weights = P7SysRDC2100$Mass)
sd.P7Sys2100 <- sqrt(var.P7Sys2100) # standard deviation of the distribution

RDCexpSysP7.2100 <- data.frame(P7Sys.2100,
                               sd.P7Sys2100,
                               row.names = "Whole System")
colnames(RDCexpSysP7.2100) <- c("Δ14C",
                                "sd")
RDCexpSysP7.2100 # summary of expected value and sd of system in a table

P7TTRDC2100 <- TTRDC(
  A = BP7,
  u = uP7,
  nyears = 1000,
  yrobs = 2100,
  C14atm = C14_Trpcs$Delta.14C,
  C14YearAD = C14_Trpcs$Year.AD,
  h = 0.1
) # for outflux, i.e. based on transit time distribution

P7TThist2100 <- C14hist(D14C = P7TTRDC2100$D14C,
                        Mass = P7TTRDC2100$Mass,
                        bin = 10) # aggregation for outflux

P7TThist2100.p <- C14hist(D14C = P7TTRDC2100$D14C,
                          Mass = P7TTRDC2100$Mass,
                          bin = 40) # with bin size b=40

P7TT.2100 <-
  weighted.mean(x = P7TTRDC2100$D14C,
                w = P7TTRDC2100$Mass) # expected value of outflux
var.P7TT2100 <-
  wtd.var(x = P7TTRDC2100$D14C, 
          weights = P7TTRDC2100$Mass)
sd.P7TT2100 <- sqrt(var.P7TT2100) # standard deviation of the distribution

RDCexpTTP7.2100 <- data.frame(P7TT.2100,
                              sd.P7TT2100,
                              row.names = "Outflux")
colnames(RDCexpTTP7.2100) <- c("Δ14C",
                               "sd")
RDCexpTTP7.2100 # summary of expected value and sd of outflux in a table

### Emanuel model
#### Age and Transit Time distributions
SAEm <- systemAge(A = BE, u = uE, a = age)
SAEmpool <- data.frame(SAEm$poolAgeDensity) # age distributions of individual pools
SAEmsys <- data.frame(SAEm$systemAgeDensity) # age distribution of the whole system
TTEm <- transitTime(A = BE, u = uE, a = age)
TTEmdens <- data.frame(TTEm$transitTimeDensity) # transit time distribution - age of outflux

#### Radiocarbon distributions
##### Year of observation = AD 1965
EmPoolRDC1965 <- PoolRDC(
  A = BE,
  u = uE,
  nyears = 1000,
  yrobs = 1965,
  C14atm = C14_NH$Delta.14C,
  C14YearAD = C14_NH$Year.AD,
  h = 0.1
) # for pools

EmP1hist1965 <- C14hist(D14C = EmPoolRDC1965$D14C,
                        Mass = EmPoolRDC1965$Mass.1,
                        bin = 40) # aggregation for pool 1
EmP2hist1965 <- C14hist(D14C = EmPoolRDC1965$D14C,
                        Mass = EmPoolRDC1965$Mass.2,
                        bin = 40) # aggregation for pool 2
EmP3hist1965 <- C14hist(D14C = EmPoolRDC1965$D14C,
                        Mass = EmPoolRDC1965$Mass.3,
                        bin = 40) # aggregation for pool 3
EmP4hist1965 <- C14hist(D14C = EmPoolRDC1965$D14C,
                        Mass = EmPoolRDC1965$Mass.4,
                        bin = 40) #aggregation for pool 4
EmP5hist1965 <- C14hist(D14C = EmPoolRDC1965$D14C,
                        Mass = EmPoolRDC1965$Mass.5,
                        bin = 40) # aggregation for pool 5

EmP1.1965 <-
  weighted.mean(x = EmPoolRDC1965$D14C,
                w = EmPoolRDC1965$Mass.1) # expected value for pool 1

EmP2.1965 <-
  weighted.mean(x = EmPoolRDC1965$D14C,
                w = EmPoolRDC1965$Mass.2) # expected value for pool 2

EmP3.1965 <-
  weighted.mean(x = EmPoolRDC1965$D14C,
                w = EmPoolRDC1965$Mass.3) # expected value for pool 3

EmP4.1965 <-
  weighted.mean(x = EmPoolRDC1965$D14C,
                w = EmPoolRDC1965$Mass.4) # expected value for pool 4

EmP5.1965 <-
  weighted.mean(x = EmPoolRDC1965$D14C,
                w = EmPoolRDC1965$Mass.5) # expected value for pool 5

RDCexpPoolsEm.1965 <- data.frame(EmP1.1965,
                                 EmP2.1965,
                                 EmP3.1965,
                                 EmP4.1965,
                                 EmP5.1965,
                                 row.names = "Δ14C")
colnames(RDCexpPoolsEm.1965) <- c(poolnamesEm[1],
                                  poolnamesEm[2],
                                  poolnamesEm[3],
                                  poolnamesEm[4],
                                  poolnamesEm[5])
RDCexpPoolsEm.1965 # summary of expected values of pools in a table

EmSysRDC1965 <- SystemRDC(
  A = BE,
  u = uE,
  nyears = 1000,
  yrobs = 1965,
  C14atm = C14_NH$Delta.14C,
  C14YearAD = C14_NH$Year.AD,
  h = 0.1
) # for the whole system

EmSyshist1965 <- C14hist(D14C = EmSysRDC1965$D14C,
                         Mass = EmSysRDC1965$Mass,
                         bin = 40) # aggregation for the whole system
EmSys.1965 <-
  weighted.mean(x = EmSysRDC1965$D14C,
                w = EmSysRDC1965$Mass) # expected value of system
var.EmSys1965 <-
  wtd.var(x = EmSysRDC1965$D14C,
          weights = EmSysRDC1965$Mass)
sd.EmSys1965 <- sqrt(var.EmSys1965) # standard deviation of the distribution

RDCexpSysEm.1965 <- data.frame(EmSys.1965,
                               sd.EmSys1965,
                               row.names = "Whole System")
colnames(RDCexpSysEm.1965) <- c("Δ14C",
                                "sd")
RDCexpSysEm.1965 # summary of expected value and sd for system in a table

EmTTRDC1965 <- TTRDC(
  A = BE,
  u = uE,
  nyears = 1000,
  yrobs = 1965,
  C14atm = C14_NH$Delta.14C,
  C14YearAD = C14_NH$Year.AD,
  h = 0.1
) # for the outflux, i.e. based on transit time distribution

EmTThist1965 <- C14hist(D14C = EmTTRDC1965$D14C,
                        Mass = EmTTRDC1965$Mass,
                        bin = 40) # aggregation for outflux

EmTT.1965 <-
  weighted.mean(x = EmTTRDC1965$D14C,
                w = EmTTRDC1965$Mass) # expected value for outflux
var.EmTT1965 <-
  wtd.var(x = EmTTRDC1965$D14C, 
          weights = EmTTRDC1965$Mass)
sd.EmTT1965 <- sqrt(var.EmTT1965) # standard deviation of the distribution

RDCexpTTEm.1965 <- data.frame(EmTT.1965,
                              sd.EmTT1965,
                              row.names = "Outflux")
colnames(RDCexpTTEm.1965) <- c("Δ14C",
                               "sd")
RDCexpTTEm.1965 # summary of expected value and sd for outflux in a table

##### Year of observation = AD 2027
EmPoolRDC2027 <- PoolRDC(
  A = BE,
  u = uE,
  nyears = 1000,
  yrobs = 2027,
  C14atm = C14_NH$Delta.14C,
  C14YearAD = C14_NH$Year.AD,
  h = 0.1
) # for pools

EmP1hist2027 <- C14hist(D14C = EmPoolRDC2027$D14C,
                        Mass = EmPoolRDC2027$Mass.1,
                        bin = 10) # aggregation for pool 1
EmP2hist2027 <- C14hist(D14C = EmPoolRDC2027$D14C,
                        Mass = EmPoolRDC2027$Mass.2,
                        bin = 10) # aggregation for pool 2
EmP3hist2027 <- C14hist(D14C = EmPoolRDC2027$D14C,
                        Mass = EmPoolRDC2027$Mass.3,
                        bin = 10) # aggregation for pool 3
EmP4hist2027 <- C14hist(D14C = EmPoolRDC2027$D14C,
                        Mass = EmPoolRDC2027$Mass.4,
                        bin = 10) # aggregation for pool 4
EmP5hist2027 <- C14hist(D14C = EmPoolRDC2027$D14C,
                        Mass = EmPoolRDC2027$Mass.5,
                        bin = 10) # aggregation for pool 5

EmP1.2027 <-
  weighted.mean(x = EmPoolRDC2027$D14C,
                w = EmPoolRDC2027$Mass.1) # expected value for pool 1
EmP2.2027 <-
  weighted.mean(x = EmPoolRDC2027$D14C,
                w = EmPoolRDC2027$Mass.2) # expected value for pool 2
EmP3.2027 <-
  weighted.mean(x = EmPoolRDC2027$D14C,
                w = EmPoolRDC2027$Mass.3) # expected value for pool 3
EmP4.2027 <-
  weighted.mean(x = EmPoolRDC2027$D14C,
                w = EmPoolRDC2027$Mass.4) # expected value for pool 4
EmP5.2027 <-
  weighted.mean(x = EmPoolRDC2027$D14C,
                w = EmPoolRDC2027$Mass.5) # expected value for pool 5

RDCexpPoolsEm.2027 <- data.frame(EmP1.2027,
                                 EmP2.2027,
                                 EmP3.2027,
                                 EmP4.2027,
                                 EmP5.2027,
                                 row.names = "Δ14C")
colnames(RDCexpPoolsEm.2027) <- c(poolnamesEm[1],
                                  poolnamesEm[2],
                                  poolnamesEm[3],
                                  poolnamesEm[4],
                                  poolnamesEm[5])
RDCexpPoolsEm.2027 # summary of expected value for pools in a table

EmSysRDC2027 <- SystemRDC(
  A = BE,
  u = uE,
  nyears = 1000,
  yrobs = 2027,
  C14atm = C14_NH$Delta.14C,
  C14YearAD = C14_NH$Year.AD,
  h = 0.1
) # for the whole system

EmSyshist2027 <- C14hist(D14C = EmSysRDC2027$D14C,
                         Mass = EmSysRDC2027$Mass,
                         bin = 10) # aggregation for system

EmSys.2027 <-
  weighted.mean(x = EmSysRDC2027$D14C,
                w = EmSysRDC2027$Mass) # expected value for system
var.EmSys2027 <-
  wtd.var(x = EmSysRDC2027$D14C,
          weights = EmSysRDC2027$Mass)
sd.EmSys2027 <- sqrt(var.EmSys2027)

RDCexpSysEm.2027 <- data.frame(EmSys.2027,
                               sd.EmSys2027,
                               row.names = "Whole System")
colnames(RDCexpSysEm.2027) <- c("Δ14C",
                                "sd")
RDCexpSysEm.2027 # summary of expected value and sd for system in a table

EmTTRDC2027 <- TTRDC(
  A = BE,
  u = uE,
  nyears = 1000,
  yrobs = 2027,
  C14atm = C14_NH$Delta.14C,
  C14YearAD = C14_NH$Year.AD,
  h = 0.1
) # for the outflux, i.e. based on transit time distribution

EmTThist2027 <- C14hist(D14C = EmTTRDC2027$D14C,
                        Mass = EmTTRDC2027$Mass,
                        bin = 10) #aggregation for outflux

EmTThist2027.p <- C14hist(D14C = EmTTRDC2027$D14C,
                          Mass = EmTTRDC2027$Mass,
                          bin = 40) # with bin size b=40

EmTT.2027 <-
  weighted.mean(x = EmTTRDC2027$D14C,
                w = EmTTRDC2027$Mass) # expected value for outflux
var.EmTT2027 <-
  wtd.var(x = EmTTRDC2027$D14C, 
          weights = EmTTRDC2027$Mass)
sd.EmTT2027 <- sqrt(var.EmTT2027) # standard deviation of the distribution

RDCexpTTEm.2027 <- data.frame(EmTT.2027,
                              sd.EmTT2027,
                              row.names = "Outflux")
colnames(RDCexpTTEm.2027) <- c("Δ14C",
                               "sd")
RDCexpTTEm.2027 # summary of expected value and sd for outflux in a table

##### Year of observation = AD 2100
EmPoolRDC2100 <- PoolRDC(
  A = BE,
  u = uE,
  nyears = 1000,
  yrobs = 2100,
  C14atm = C14_NH$Delta.14C,
  C14YearAD = C14_NH$Year.AD,
  h = 0.1
) # for pools

EmP1hist2100 <- C14hist(D14C = EmPoolRDC2100$D14C,
                        Mass = EmPoolRDC2100$Mass.1,
                        bin = 10) # aggregation for pool 1
EmP2hist2100 <- C14hist(D14C = EmPoolRDC2100$D14C,
                        Mass = EmPoolRDC2100$Mass.2,
                        bin = 10) # aggregation for pool 2
EmP3hist2100 <- C14hist(D14C = EmPoolRDC2100$D14C,
                        Mass = EmPoolRDC2100$Mass.3,
                        bin = 10) # aggregation for pool 3
EmP4hist2100 <- C14hist(D14C = EmPoolRDC2100$D14C,
                        Mass = EmPoolRDC2100$Mass.4,
                        bin = 10) # aggregation for pool 4
EmP5hist2100 <- C14hist(D14C = EmPoolRDC2100$D14C,
                        Mass = EmPoolRDC2100$Mass.5,
                        bin = 10) # aggregation for pool 5

EmP1.2100 <-
  weighted.mean(x = EmPoolRDC2100$D14C,
                w = EmPoolRDC2100$Mass.1) # expected value for pool 1
EmP2.2100 <-
  weighted.mean(x = EmPoolRDC2100$D14C,
                w = EmPoolRDC2100$Mass.2) # expected value for pool 2
EmP3.2100 <-
  weighted.mean(x = EmPoolRDC2100$D14C,
                w = EmPoolRDC2100$Mass.3) # expected value for pool 3
EmP4.2100 <-
  weighted.mean(x = EmPoolRDC2100$D14C,
                w = EmPoolRDC2100$Mass.4) # expected value for pool 4
EmP5.2100 <-
  weighted.mean(x = EmPoolRDC2100$D14C,
                w = EmPoolRDC2100$Mass.5) # expected value for pool 5

RDCexpPoolsEm.2100 <- data.frame(EmP1.2100,
                                 EmP2.2100,
                                 EmP3.2100,
                                 EmP4.2100,
                                 EmP5.2100,
                                 row.names = "Δ14C")
colnames(RDCexpPoolsEm.2100) <- c(poolnamesEm[1],
                                  poolnamesEm[2],
                                  poolnamesEm[3],
                                  poolnamesEm[4],
                                  poolnamesEm[5])
RDCexpPoolsEm.2100 # summary of expected values for pools in a table

EmSysRDC2100 <- SystemRDC(
  A = BE,
  u = uE,
  nyears = 1000,
  yrobs = 2100,
  C14atm = C14_NH$Delta.14C,
  C14YearAD = C14_NH$Year.AD,
  h = 0.1
) # for the whole system

EmSyshist2100 <- C14hist(D14C = EmSysRDC2100$D14C,
                         Mass = EmSysRDC2100$Mass,
                         bin = 10) # aggregation for system

EmSys.2100 <-
  weighted.mean(x = EmSysRDC2100$D14C,
                w = EmSysRDC2100$Mass) # expected value for system
var.EmSys2100 <-
  wtd.var(x = EmSysRDC2100$D14C,
          weights = EmSysRDC2100$Mass)
sd.EmSys2100 <- sqrt(var.EmSys2100) # standard deviation of the distribution

RDCexpSysEm.2100 <- data.frame(EmSys.2100,
                               sd.EmSys2100,
                               row.names = "Whole System")
colnames(RDCexpSysEm.2100) <- c("Δ14C",
                                "sd")
RDCexpSysEm.2100 # summary of expected value and sd of system in a table

EmTTRDC2100 <- TTRDC(
  A = BE,
  u = uE,
  nyears = 1000,
  yrobs = 2100,
  C14atm = C14_NH$Delta.14C,
  C14YearAD = C14_NH$Year.AD,
  h = 0.1
) # for the outflux, i.e. based on transit time distribution

EmTThist2100 <- C14hist(D14C = EmTTRDC2100$D14C,
                        Mass = EmTTRDC2100$Mass,
                        bin = 10) # aggregation for outflux

EmTThist2100.p <- C14hist(D14C = EmTTRDC2100$D14C,
                          Mass = EmTTRDC2100$Mass,
                          bin = 40) # with bin size b=40

EmTT.2100 <-
  weighted.mean(x = EmTTRDC2100$D14C,
                w = EmTTRDC2100$Mass) # expected value of outflux
var.EmTT2100 <-
  wtd.var(x = EmTTRDC2100$D14C, 
          weights = EmTTRDC2100$Mass)
sd.EmTT2100 <- sqrt(var.EmTT2100) # standard deviation of the distribution

RDCexpTTEm.2100 <- data.frame(EmTT.2100,
                              sd.EmTT2100,
                              row.names = "Outflux")
colnames(RDCexpTTEm.2100) <- c("Δ14C",
                               "sd")
RDCexpTTEm.2100 # summary of expected value and sd of outflux in a table

## Plots
### Figure 3 - Scheme of merge of 14C curves for the period -53050 to 2100 and data sets used

nhtrpcs.curves <- data.frame(C14_Trpcs, C14_NH) # column names: none=Tropics, 1=NH

# C14nhtrpcs.p <- 
ggplot(data = nhtrpcs.curves) +
  geom_point(
    mapping = aes(x = Year.AD,
                  y = Delta.14C,
                  color = "Tropics"),
    alpha = .6,
    size = 2
  ) +
  geom_point(
    mapping = aes(x = Year.AD.1,
                  y = Delta.14C.1,
                  color = "C14NH"),
    alpha = .6,
    size = 2
  ) +
  labs(x = "Calendar Years (AD)",
       y = expression(paste(Delta ^ {
         14
       }, "C (\u2030)"))) +
  scale_color_manual(
    name = "Radiocarbon curves",
    values = c("Tropics" = "darkolivegreen4",
               "C14NH" = "deeppink3"),
    labels = c("Tropics" = "Tropics",
               "C14NH" = "Northern Hemisphere")
  ) +
  ggtitle("Radiocarbon curves") +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(face = "bold")
  ) +
  coord_cartesian(xlim = c(1850, 2100)) +
  annotate(
    "rect",
    xmin = 1830,
    xmax = 1950,
    ymin = -300,
    ymax = 860,
    alpha = .4,
    fill = "grey",
    color = "black"
  ) +
  annotate("label",
           x = 1895,
           y = -150,
           label = "IntCal20\n(Reimer et al., 2020)") +
  annotate(
    "rect",
    xmin = 1950,
    xmax = 2015,
    ymin = -300,
    ymax = 860,
    alpha = .1,
    fill = "grey",
    color = "black"
  ) +
  annotate("label",
           x = 1983,
           y = -150,
           label = "Records\n(Graven et al., 2017)"
  ) +
  annotate(
    "rect",
    xmin = 2015,
    xmax = 2100,
    ymin = -300,
    ymax = 860,
    alpha = .4,
    fill = "grey",
    color = "black"
  ) +
  annotate("label",
           x = 2060,
           y = -150,
           label = "Forecast RCP8.5\n(Graven, 2015)")

### Figure 4
# pools' age distributions
poolHages <-
  ggplot(data = SAHpool,
         mapping = aes(x = age)) +
  geom_line(mapping = aes(y = X1, 
                          color = "P1"),
            size = 1.1) +
  geom_line(mapping = aes(y = X2, 
                          color = "P2"),
            size = 1.1) +
  geom_line(mapping = aes(y = X3,
                          color = "P3"),
            size = 1.1) +
  geom_line(mapping = aes(y = X4, 
                          color = "P4"),
            size = 1.1) +
  geom_line(mapping = aes(y = X5,
                          color = "P5"),
            size = 1.1) +
  geom_line(mapping = aes(y = X6,
                          color = "P6"),
            size = 1.1) +
  geom_line(mapping = aes(y = X7,
                          color = "P7"),
            size = 1.1) +
  coord_cartesian(xlim = c(0, 150),
                  ylim = c(0, 0.2)) +
  labs(x = "Age (years)",
       y = "Probability density",
       color = "Pools") +
  scale_color_manual(
    labels = c(
      "P1" = bquote(paste(x[1], ": Dead Roots")),
      "P2" = bquote(paste(x[2], ": Oi")),
      "P3" = bquote(paste(x[3], ": Oe/a L")),
      "P4" = bquote(paste(x[4], ": Oe/a H")),
      "P5" = bquote(paste(x[5], ": A, LF (> 80 µm)")),
      "P6" = bquote(paste(x[6], ": A, LF (< 80 µm)")),
      "P7" = bquote(paste(x[7], ": Mineral Associated"))
    ),
    values = c(
      "P1" = "darkolivegreen3",
      "P2" = "blue",
      "P3" = "#009E73",
      "P4" = "#F0E442",
      "P5" = "deeppink3",
      "P6" = "#FE9929",
      "P7" = "#CC79A7"
    )
  ) +
  theme(legend.position = c(0.75, 0.65))

# AD 1965 - pools' 14C distributions - stacked
PoolH_RDC1965 <- ggplot(mapping = aes(x = mid,
                                      y = M)) +
  geom_bar(
    data = HP1hist1965,
    mapping = aes(fill = "P1"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = HP2hist1965,
    mapping = aes(fill = "P2"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = HP3hist1965,
    mapping = aes(fill = "P3"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = HP4hist1965,
    mapping = aes(fill = "P4"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = HP5hist1965,
    mapping = aes(fill = "P5"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = HP6hist1965,
    mapping = aes(fill = "P6"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = HP7hist1965,
    mapping = aes(fill = "P7"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  )

LegPoolH1965 <- PoolH_RDC1965 +
  labs(y = expression(paste("Mass of Carbon (g.",m^-2,")")),
       x = expression(paste(Delta ^ {14}, "C (\u2030)")),
       fill = "Pools") +
  scale_fill_manual(
    labels = c(
      "P1" = poolnamesHFS[1],
      "P2" = poolnamesHFS[2],
      "P3" = poolnamesHFS[3],
      "P4" = poolnamesHFS[4],
      "P5" = poolnamesHFS[5],
      "P6" = poolnamesHFS[6],
      "P7" = poolnamesHFS[7]
    ),
    values = c(
      "P1" = "darkolivegreen3",
      "P2" = "blue",
      "P3" = "#009E73",
      "P4" = "#F0E442",
      "P5" = "deeppink3",
      "P6" = "#FE9929",
      "P7" = "#CC79A7"
    )
  ) +
  ggtitle("Year of observation = 1965") +
  theme(legend.position = c(0.75, 0.65)) +
  coord_cartesian(xlim = c(-300, 800))

# AD 2027 - pools' 14C distributions - stacked
PoolH_RDC2027 <- ggplot(mapping = aes(x = mid,
                                      y = M)) +
  geom_bar(
    data = HP1hist2027,
    mapping = aes(fill = "P1"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = HP2hist2027,
    mapping = aes(fill = "P2"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = HP3hist2027,
    mapping = aes(fill = "P3"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = HP4hist2027,
    mapping = aes(fill = "P4"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = HP5hist2027,
    mapping = aes(fill = "P5"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = HP6hist2027,
    mapping = aes(fill = "P6"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = HP7hist2027,
    mapping = aes(fill = "P7"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  )

LegPoolH2027 <- PoolH_RDC2027 +
  labs(y = expression(paste("Mass of Carbon (g.",m^-2,")")),
       x = expression(paste(Delta ^ {14}, "C (\u2030)")),
       fill = "Pools") +
  scale_fill_manual(
    labels = c(
      "P1" = poolnamesHFS[1],
      "P2" = poolnamesHFS[2],
      "P3" = poolnamesHFS[3],
      "P4" = poolnamesHFS[4],
      "P5" = poolnamesHFS[5],
      "P6" = poolnamesHFS[6],
      "P7" = poolnamesHFS[7]
    ),
    values = c(
      "P1" = "darkolivegreen3",
      "P2" = "blue",
      "P3" = "#009E73",
      "P4" = "#F0E442",
      "P5" = "deeppink3",
      "P6" = "#FE9929",
      "P7" = "#CC79A7"
    )
  ) +
  ggtitle("Year of observation = 2027") +
  theme(legend.position = c(0.75, 0.65)) +
  coord_cartesian(xlim = c(-300, 800))

# AD 2100 - pools' 14C distributions - stacked
PoolH_RDC2100 <- ggplot(mapping = aes(x = mid,
                                      y = M)) +
  geom_bar(
    data = HP1hist2100,
    mapping = aes(fill = "P1"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = HP2hist2100,
    mapping = aes(fill = "P2"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = HP3hist2100,
    mapping = aes(fill = "P3"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = HP4hist2100,
    mapping = aes(fill = "P4"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = HP5hist2100,
    mapping = aes(fill = "P5"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = HP6hist2100,
    mapping = aes(fill = "P6"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = HP7hist2100,
    mapping = aes(fill = "P7"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  )

LegPoolH2100 <- PoolH_RDC2100 +
  labs(y = expression(paste("Mass of Carbon (g.",m^-2,")")),
       x = expression(paste(Delta ^ {14}, "C (\u2030)")),
       fill = "Pools") +
  scale_fill_manual(
    labels = c(
      "P1" = poolnamesHFS[1],
      "P2" = poolnamesHFS[2],
      "P3" = poolnamesHFS[3],
      "P4" = poolnamesHFS[4],
      "P5" = poolnamesHFS[5],
      "P6" = poolnamesHFS[6],
      "P7" = poolnamesHFS[7]
    ),
    values = c(
      "P1" = "darkolivegreen3",
      "P2" = "blue",
      "P3" = "#009E73",
      "P4" = "#F0E442",
      "P5" = "deeppink3",
      "P6" = "#FE9929",
      "P7" = "#CC79A7"
    )
  ) +
  ggtitle("Year of observation = 2100") +
  theme(legend.position = c(0.75, 0.65)) +
  coord_cartesian(xlim = c(-300, 800))

# arranging alltogether
Hallp <- ggpubr::ggarrange(
  poolHages +
    ggtitle("Any year of observation") +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 24,
                                    face = "bold"),
          axis.title = element_text(size = 28),
          axis.text = element_text(size = 28),
          legend.text = element_text(size = 24),
          legend.title = element_text(size = 24)),
  LegPoolH1965 +
    ggtitle("Year of observation = 1965") +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 24,
                                    face = "bold"),
          axis.title = element_text(size = 28),
          axis.text = element_text(size = 28),
          legend.text = element_text(size = 24),
          legend.title = element_text(size = 24)),
  LegPoolH2027 +
    ggtitle("Year of observation = 2027") +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 24,
                                    face = "bold"),
          axis.title = element_text(size = 28),
          axis.text = element_text(size = 28),
          legend.text = element_text(size = 24),
          legend.title = element_text(size = 24)),
  LegPoolH2100 +
    ggtitle("Year of observation = 2100") +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 24,
                                    face = "bold"),
          axis.title = element_text(size = 28),
          axis.text = element_text(size = 28),
          legend.text = element_text(size = 24),
          legend.title = element_text(size = 24)),
  # list of plots
  labels = "auto",
  # labels
  common.legend = T,
  # common legend
  legend = "right",
  # legend position
  align = "hv",
  # Align them both, horizontal and vertical
  nrow = 2, # number of rows
  ncol = 2,  # number of columns
  font.label = list(size = 24)) 

H.allplot <- annotate_figure(Hallp,
                             top = text_grob(
                               "HARVARD FOREST SOIL MODEL",
                               color = "black",
                               size = 28
                             ))
# Final - Figure 4 - Distributions of age for each pool in HFS model and
# 14C distributions of pools in 1965, 2027 and 2100
H.allplot

### Figure 5
# transit time distribution
TTHdist <- ggplot(data = TTHdens,
                  mapping = aes(x = age)) +
  geom_line(mapping = aes(y = TTH.transitTimeDensity, 
                          color = "TT"),
            size = 1.3) +
  coord_cartesian(xlim = c(0, 150),
                  ylim = c(0, 0.1)) +
  labs(x = "Age (years)",
       y = "Probability density",
       color = element_blank()) +
  scale_color_manual(labels = c("TT" = "Outflux"),
                     values = c("TT" = "#F8766D")) +
  theme(legend.position = c(0.75, 0.65))

# 14C distributions
allTTHFSRDC <-
  ggplot(mapping = aes(x = mid,
                       y = M)) +
  geom_bar(
    data = HTThist1965,
    mapping = aes(color = "y1",
                  fill = "y1"),
    stat = "identity",
    position = position_dodge(),
    alpha = .4
  ) +
  geom_bar(
    data = HTThist2027.p,
    mapping = aes(color = "y2",
                  fill = "y2"),
    stat = "identity",
    position = position_dodge(),
    alpha = .3
  ) +
  geom_bar(
    data = HTThist2100.p,
    mapping = aes(color = "y3",
                  fill = "y3"),
    stat = "identity",
    position = position_dodge(),
    alpha = .3
  ) +
  labs(y = expression(paste("Mass of Carbon (g.", m ^ -2, ")")),
       x = expression(paste(Delta ^ {
         14
       }, "C (\u2030)"))) +
  scale_color_manual(
    name = "Year of Observation",
    labels = c("y1" = "1965 AD",
               "y2" = "2027 AD",
               "y3" = "2100 AD"),
    values = c("y1" = "blue",
               "y2" = "darkolivegreen3",
               "y3" = "deeppink3")
  )  +
  scale_fill_manual(
    name = "Year of Observation",
    labels = c("y1" = "1965 AD",
               "y2" = "2027 AD",
               "y3" = "2100 AD"),
    values = c("y1" = "blue",
               "y2" = "darkolivegreen3",
               "y3" = "deeppink3")
  ) +
  theme(legend.position = c(0.8, 0.65))

HTTallp <- ggpubr::ggarrange(
  TTHdist +
    ggtitle("Transit time distribution", subtitle = "Any year of observation") +
    theme(
      plot.title = element_text(hjust = 0.5,
                                size = 24,
                                face = "bold"),
      legend.position = "none",
      axis.title = element_text(size = 28),
      axis.text = element_text(size = 28),
      legend.text = element_text(size = 24),
      legend.title = element_text(size = 24),
      plot.subtitle = element_text(size = 24)
    ),
  allTTHFSRDC +
    ggtitle("Outflux radiocarbon distribution") +
    theme(
      plot.title = element_text(hjust = 0.5,
                                size = 24,
                                face = "bold"),
      legend.position = "top",
      axis.title = element_text(size = 28),
      axis.text = element_text(size = 28),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16)
    ),
  # list of plots
  labels = "auto",
  # labels
  common.legend = F,
  # common legend
  align = "v",
  # Align them both, horizontal and vertical
  nrow = 1,
  # number of rows
  font.label = list(size = 24),
  ncol = 2  # number of columns
) 

HTT.all <- annotate_figure(HTTallp,
                           top = text_grob(
                             "HARVARD FOREST SOIL MODEL",
                             color = "black",
                             size = 28))
# Final Figure 5 - HFS transit time distribution and 14C distributions
# for the years 1965, 2027 and 2100
HTT.all

### Figure 6 
resp_yr <- HSMdf %>% 
  mutate(Year = floor(Year)) %>% 
  arrange(desc(Year))

resp1996 <- resp_yr %>% 
  filter(Year == 1996)
resp1998 <- resp_yr %>% 
  filter(Year == 1998)
resp2002 <- resp_yr %>% 
  filter(Year == 2002)
resp2008 <- resp_yr %>% 
  filter(Year == 2008)

# plots from measurements
# AD 1996
exp1996 <- ggplot(data = resp1996, aes(x = D14C)) +
  geom_histogram(binwidth = 10,
                 fill = "grey66",
                 color = "black") +
  geom_vline(
    aes(xintercept = mean(D14C), color = "mean1996"),
    linetype = "dashed",
    size = 1.5
  ) +
  labs(y = "Frequency",
       x = expression(paste(Delta ^ {
         14
       }, "C (\u2030)"))) +
  scale_color_manual(name = element_blank(),
                     values = c("mean1996" = "#D55E00"),
                     labels = c("mean1996" = expression(paste("Mean ", Delta^{14}, "C = 129.6 \u2030")))) +
  coord_cartesian(x = c(-50, 500)) +
  ggtitle(label = expression(paste(Delta ^ {
    14
  }, "C in total soil C", O[2], " efflux")), subtitle = "Measurements from 1996") +
  theme(legend.position = c(0.7, 0.6)) +
  geom_text(label = "n = 12",
            x = 0,
            y = 1.5,
            size = 4) # to save plot use size = 7
# AD 1998
exp1998 <- ggplot(data = resp1998, aes(x = D14C)) +
  geom_histogram(binwidth = 10,
                 fill = "grey66",
                 color = "black") +
  geom_vline(
    aes(xintercept = mean(D14C), color = "mean1998"),
    linetype = "dashed",
    size = 1.5
  ) +
  labs(y = "Frequency",
       x = expression(paste(Delta ^ {
         14
       }, "C (\u2030)"))) +
  scale_color_manual(name = expression(paste(Delta ^ {
    14
  }, "C (\u2030)")),
  values = c("mean1998" = "#D55E00"),
  labels = c("mean1998" = expression(paste("Mean ", Delta^{14}, "C = 117.6 \u2030")))) +
  coord_cartesian(x = c(-50, 500)) +
  ggtitle(label = expression(paste(Delta ^ {
    14
  }, "C in total soil C", O[2], " efflux")), subtitle = "Measurements from 1998") +
  theme(
    legend.position = c(0.85, 0.6)) +
  geom_text(label = "n = 28",
            x = 0,
            y = 3,
            size = 4) # to save plot use size = 7
# AD 2002
exp2002 <- ggplot(data = resp2002, aes(x = D14C)) +
  geom_histogram(binwidth = 10,
                 fill = "grey66",
                 color = "black") +
  geom_vline(
    aes(xintercept = mean(D14C), color = "mean2002"),
    linetype = "dashed",
    size = 1.5
  ) +
  labs(y = "Frequency",
       x = expression(paste(Delta ^ {
         14
       }, "C (\u2030)"))) +
  scale_color_manual(name = expression(paste(Delta ^ {
    14
  }, "C (\u2030)")),
  values = c("mean2002" = "#D55E00"),
  labels = c("mean2002" = expression(paste("Mean ", Delta^{14}, "C = 100.8 \u2030")))) +
  coord_cartesian(x = c(-50, 500)) +
  ggtitle(label = expression(paste(Delta ^ {
    14
  }, "C in total soil C", O[2], " efflux")), subtitle = "Measurements from 2002") +
  theme(legend.position = c(0.7, 0.6)) +
  geom_text(label = "n = 23",
            x = 0,
            y = 3,
            size = 4) # to save plot use size = 7
# AD 2008
exp2008 <- ggplot(data = resp2008, aes(x = D14C)) +
  geom_histogram(binwidth = 10,
                 fill = "grey66",
                 color = "black") +
  geom_vline(
    aes(xintercept = mean(D14C), color = "mean2008"),
    linetype = "dashed",
    size = 1.5
  ) +
  labs(y = "Frequency",
       x = expression(paste(Delta ^ {
         14
       }, "C (\u2030)"))) +
  scale_color_manual(name = expression(paste(Delta ^ {
    14
  }, "C (\u2030)")),
  values = c("mean2008" = "#D55E00"),
  labels = c("mean2008" = expression(paste("Mean ", Delta^{14}, "C = 74.9 \u2030")))) +
  coord_cartesian(x = c(-50, 500)) +
  ggtitle(label = expression(paste(Delta ^ {
    14
  }, "C in total soil C", O[2], " efflux")), subtitle = "Measurements from 2008") +
  theme(legend.position = c(0.7, 0.6)) +
  geom_text(label = "n = 10",
            x = 0,
            y = 1,
            size = 4) # to save plot use size = 7
# from theoretical estimations
# AD 1996
TTRDC1996.H <- ggplot(data = HTThist1996,
                      mapping = aes(x = mid,
                                    y = M)) +
  geom_bar(stat = "identity",
           position = position_dodge(),
           fill = "grey88",
           color = "black") +
  geom_vline(aes(xintercept = HTT.1996,
                 color = "HTT.1996"),
             linetype = "dashed") + 
  labs(y = expression(paste("Mass of Carbon (g.",m^-2,")")),
       x = expression(paste(Delta ^ {14}, "C (\u2030)"))) +
  scale_color_manual(name = element_blank(),
                     labels = c("HTT.1996" = expression(paste("Expected ", Delta^{14}, "C = 153 \u2030"))),
                     values = c("HTT.1996" = "#0072B2")) +
  ggtitle("Outflux", subtitle = "Year of observation = 1996") +
  theme(legend.position = c(0.7, 0.6)) +
  coord_cartesian(xlim = c(-50, 500))
# AD 1998
TTRDC1998.H <- ggplot(data = HTThist1998,
                      mapping = aes(x = mid,
                                    y = M)) +
  geom_bar(stat = "identity",
           position = position_dodge(),
           fill = "grey88",
           color = "black") +
  geom_vline(aes(xintercept = HTT.1998,
                 color = "HTT.1998"),
             linetype = "dashed") + 
  labs(y = expression(paste("Mass of Carbon (g.",m^-2,")")),
       x = expression(paste(Delta ^ {14}, "C (\u2030)"))) +
  scale_color_manual(name = element_blank(),
                     labels = c("HTT.1998" = expression(paste("Expected ", Delta^{14}, "C = 139.4 \u2030"))),
                     values = c("HTT.1998" = "#0072B2")) +
  ggtitle("Outflux", subtitle = "Year of observation = 1998") +
  theme(legend.position = c(0.7, 0.6)) +
  coord_cartesian(xlim = c(-50, 500))
# AD 2002
TTRDC2002.H <- ggplot(data = HTThist2002,
                      mapping = aes(x = mid,
                                    y = M)) +
  geom_bar(stat = "identity",
           position = position_dodge(),
           fill = "grey88",
           color = "black") +
  geom_vline(aes(xintercept = HTT.2002,
                 color = "HTT.2002"),
             linetype = "dashed") + 
  labs(y = expression(paste("Mass of Carbon (g.",m^-2,")")),
       x = expression(paste(Delta ^ {14}, "C (\u2030)"))) +
  scale_color_manual(name = element_blank(),
                     labels = c("HTT.2002" = expression(paste("Expected ", Delta^{14}, "C = 115.9 \u2030"))),
                     values = c("HTT.2002" = "#0072B2")) +
  ggtitle("Outflux", subtitle = "Year of observation = 2002") +
  theme(legend.position = c(0.7, 0.6)) +
  coord_cartesian(xlim = c(-50, 500))
# AD 2008
TTRDC2008.H <- ggplot(data = HTThist2008,
                      mapping = aes(x = mid,
                                    y = M)) +
  geom_bar(stat = "identity",
           position = position_dodge(),
           fill = "grey88",
           color = "black") +
  geom_vline(aes(xintercept = HTT.2008,
                 color = "HTT.2008"),
             linetype = "dashed") + 
  labs(y = expression(paste("Mass of Carbon (g.",m^-2,")")),
       x = expression(paste(Delta ^ {14}, "C (\u2030)"))) +
  scale_color_manual(name = element_blank(),
                     labels = c("HTT.2008" = expression(paste("Expected ", Delta^{14}, "C = 85 \u2030"))),
                     values = c("HTT.2008" = "#0072B2")) +
  ggtitle("Outflux", subtitle = "Year of observation = 2008") +
  theme(legend.position = c(0.7, 0.6)) +
  coord_cartesian(xlim = c(-50, 500))
# alltogether - commented sizes are the ones used for saving fig
# with ggsave: device = cairo_pdf, height = 20, width = 22, dpi = 600, units = "in"
tt_exp1996 <-
  cowplot::plot_grid(TTRDC1996.H 
                     # +
                     #   theme(
                     #     title = element_text(size = 24),
                     #     axis.title.y = element_text(size = 26),
                     #     legend.text = element_text(size = 26),
                     #     axis.text = element_text(size = 24)
                     #   )
                     , 
                     exp1996 
                     # +
                     #   theme(
                     #     title = element_text(size = 24),
                     #     axis.title.y = element_text(size = 26),
                     #     legend.text = element_text(size = 26),
                     #     axis.text = element_text(size = 24)
                     #   )
                     , 
                     align = "v",
                     ncol = 1)
tt_exp1998 <-
  cowplot::plot_grid(TTRDC1998.H 
                     # +
                     #   theme(title = element_text(size = 24),
                     #     axis.title.y = element_text(size = 26),
                     #     legend.text = element_text(size = 26),
                     #     axis.text = element_text(size = 24))
                     , 
                     exp1998 
                     # +
                     #   theme(
                     #     title = element_text(size = 24),
                     #     axis.title.y = element_text(size = 26),
                     #     legend.text = element_text(size = 26),
                     #     axis.text = element_text(size = 24)
                     #   )
                     , 
                     align = "v",
                     ncol = 1)
tt_exp2002 <-
  cowplot::plot_grid(TTRDC2002.H 
                     # +
                     #   theme(
                     #     title = element_text(size = 24),
                     #     axis.title.y = element_text(size = 26),
                     #     legend.text = element_text(size = 26),
                     #     axis.text = element_text(size = 24)
                     #   )
                     , 
                     exp2002 
                     # +
                     #   theme(
                     #     title = element_text(size = 24),
                     #     axis.title.y = element_text(size = 26),
                     #     legend.text = element_text(size = 26),
                     #     axis.text = element_text(size = 24)
                     #   )
                     , 
                     align = "v",
                     ncol = 1)
tt_exp2008 <-
  cowplot::plot_grid(TTRDC2008.H 
                     # +
                     #   theme(
                     #     title = element_text(size = 24),
                     #     axis.title.y = element_text(size = 26),
                     #     legend.text = element_text(size = 26),
                     #     axis.text = element_text(size = 24)
                     #   )
                     , 
                     exp2008 
                     # +
                     #   theme(
                     #     title = element_text(size = 24),
                     #     axis.title.y = element_text(size = 26),
                     #     legend.text = element_text(size = 26),
                     #     axis.text = element_text(size = 24)
                     #   )
                     , 
                     align = "v",
                     ncol = 1)
# Final Figure 6 - Comparison theoretical distributions of Delta14C and
# soil Delta14CO2 efflux data
cowplot::plot_grid(tt_exp1996,
                   tt_exp1998,
                   tt_exp2002,
                   tt_exp2008, 
                   align = "v",
                   ncol = 2,
                   nrow = 2,
                   labels = c("a",
                              "b",
                              "c",
                              "d"),
                   label_size = 18) # to save plot use label_size = 28

### Figure A2
# age distributions
poolP7ages <-
  ggplot(data = SAP7pool,
         mapping = aes(x = age)) +
  geom_line(mapping = aes(y = X1, 
                          color = "P1"),
            size = 1.1) +
  geom_line(mapping = aes(y = X2, 
                          color = "P2"),
            size = 1.1) +
  geom_line(mapping = aes(y = X3,
                          color = "P3"),
            size = 1.1) +
  geom_line(mapping = aes(y = X4, 
                          color = "P4"),
            size = 1.1) +
  geom_line(mapping = aes(y = X5,
                          color = "P5"),
            size = 1.1) +
  geom_line(mapping = aes(y = X6,
                          color = "P6"),
            size = 1.1) +
  geom_line(mapping = aes(y = X7,
                          color = "P7"),
            size = 1.1) +
  coord_cartesian(xlim = c(0, 25),
                  ylim = c(0, 1)) +
  labs(x = "Age (years)",
       y = "Probability density",
       color = "Pools") +
  scale_color_manual(
    labels = c(
      "P1" = bquote(paste(x[1], ": Foliage")),
      "P2" = bquote(paste(x[2], ": Wood")),
      "P3" = bquote(paste(x[3], ": Fine roots")),
      "P4" = bquote(paste(x[4], ": Coarse roots")),
      "P5" = bquote(paste(x[5], ": Fine litter")),
      "P6" = bquote(paste(x[6], ": Coarse woody debris")),
      "P7" = bquote(paste(x[7], ": Soil carbon (0 - 30 cm)"))
    ),
    values = c(
      "P1" = "darkolivegreen3",
      "P2" = "blue",
      "P3" = "#009E73",
      "P4" = "#F0E442",
      "P5" = "deeppink3",
      "P6" = "#FE9929",
      "P7" = "#CC79A7"
    )
  ) +
  theme(legend.position = c(0.75, 0.65))
# pools' 14C distributions
# AD 1965
PoolP7_RDC1965 <- ggplot(mapping = aes(x = mid,
                                       y = M)) +
  geom_bar(
    data = P7P1hist1965,
    mapping = aes(fill = "P1"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = P7P2hist1965,
    mapping = aes(fill = "P2"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = P7P3hist1965,
    mapping = aes(fill = "P3"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = P7P4hist1965,
    mapping = aes(fill = "P4"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = P7P5hist1965,
    mapping = aes(fill = "P5"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = P7P6hist1965,
    mapping = aes(fill = "P6"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = P7P7hist1965,
    mapping = aes(fill = "P7"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  )

LegPoolP71965 <- PoolP7_RDC1965 +
  labs(y = expression(paste("Mass of Carbon (Mg h",a^-1,")")),
       x = expression(paste(Delta ^ {14}, "C (\u2030)")),
       fill = "Pools") +
  scale_fill_manual(
    labels = c(
      "P1" = poolnamesP7[1],
      "P2" = poolnamesP7[2],
      "P3" = poolnamesP7[3],
      "P4" = poolnamesP7[4],
      "P5" = poolnamesP7[5],
      "P6" = poolnamesP7[6],
      "P7" = poolnamesP7[7]
    ),
    values = c(
      "P1" = "darkolivegreen3",
      "P2" = "blue",
      "P3" = "#009E73",
      "P4" = "#F0E442",
      "P5" = "deeppink3",
      "P6" = "#FE9929",
      "P7" = "#CC79A7"
    )
  ) +
  theme(legend.position = c(0.75, 0.65))
# AD 2027
PoolP7_RDC2027 <- ggplot(mapping = aes(x = mid,
                                       y = M)) +
  geom_bar(
    data = P7P1hist2027,
    mapping = aes(fill = "P1"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = P7P2hist2027,
    mapping = aes(fill = "P2"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = P7P3hist2027,
    mapping = aes(fill = "P3"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = P7P4hist2027,
    mapping = aes(fill = "P4"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = P7P5hist2027,
    mapping = aes(fill = "P5"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = P7P6hist2027,
    mapping = aes(fill = "P6"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = P7P7hist2027,
    mapping = aes(fill = "P7"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  )

LegPoolP72027 <- PoolP7_RDC2027 +
  labs(y = expression(paste("Mass of Carbon (Mg h",a^-1,")")),
       x = expression(paste(Delta ^ {14}, "C (\u2030)")),
       fill = "Pools") +
  scale_fill_manual(
    labels = c(
      "P1" = poolnamesP7[1],
      "P2" = poolnamesP7[2],
      "P3" = poolnamesP7[3],
      "P4" = poolnamesP7[4],
      "P5" = poolnamesP7[5],
      "P6" = poolnamesP7[6],
      "P7" = poolnamesP7[7]
    ),
    values = c(
      "P1" = "darkolivegreen3",
      "P2" = "blue",
      "P3" = "#009E73",
      "P4" = "#F0E442",
      "P5" = "deeppink3",
      "P6" = "#FE9929",
      "P7" = "#CC79A7"
    )
  ) +
  theme(legend.position = c(0.75, 0.65))
#AD 2100
PoolP7_RDC2100 <- ggplot(mapping = aes(x = mid,
                                       y = M)) +
  geom_bar(
    data = P7P1hist2100,
    mapping = aes(fill = "P1"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = P7P2hist2100,
    mapping = aes(fill = "P2"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = P7P3hist2100,
    mapping = aes(fill = "P3"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = P7P4hist2100,
    mapping = aes(fill = "P4"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = P7P5hist2100,
    mapping = aes(fill = "P5"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = P7P6hist2100,
    mapping = aes(fill = "P6"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = P7P7hist2100,
    mapping = aes(fill = "P7"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  )

LegPoolP72100 <- PoolP7_RDC2100 +
  labs(y = expression(paste("Mass of Carbon (Mg h",a^-1,")")),
       x = expression(paste(Delta ^ {14}, "C (\u2030)")),
       fill = "Pools") +
  scale_fill_manual(
    labels = c(
      "P1" = poolnamesP7[1],
      "P2" = poolnamesP7[2],
      "P3" = poolnamesP7[3],
      "P4" = poolnamesP7[4],
      "P5" = poolnamesP7[5],
      "P6" = poolnamesP7[6],
      "P7" = poolnamesP7[7]
    ),
    values = c(
      "P1" = "darkolivegreen3",
      "P2" = "blue",
      "P3" = "#009E73",
      "P4" = "#F0E442",
      "P5" = "deeppink3",
      "P6" = "#FE9929",
      "P7" = "#CC79A7"
    )
  ) +
  theme(legend.position = c(0.75, 0.65))

# arranging alltogether
P7allp <- ggpubr::ggarrange(
  poolP7ages +
    ggtitle("Any year of observation") +
    theme(
      plot.title = element_text(hjust = 0.5,
                                size = 24,
                                face = "bold"),
      axis.title = element_text(size = 28),
      axis.text = element_text(size = 28),
      legend.text = element_text(size = 24),
      legend.title = element_text(size = 24)
    ),
  LegPoolP71965 +
    ggtitle("Year of observation = 1965") +
    theme(
      plot.title = element_text(hjust = 0.5,
                                size = 24,
                                face = "bold"),
      axis.title = element_text(size = 28),
      axis.text = element_text(size = 28),
      legend.text = element_text(size = 24),
      legend.title = element_text(size = 24)
    ),
  LegPoolP72027 +
    ggtitle("Year of observation = 2027") +
    theme(
      plot.title = element_text(hjust = 0.5,
                                size = 24,
                                face = "bold"),
      axis.title = element_text(size = 28),
      axis.text = element_text(size = 28),
      legend.text = element_text(size = 24),
      legend.title = element_text(size = 24)
    ),
  LegPoolP72100 +
    ggtitle("Year of observation = 2100") +
    theme(
      plot.title = element_text(hjust = 0.5,
                                size = 24,
                                face = "bold"),
      axis.title = element_text(size = 28),
      axis.text = element_text(size = 28),
      legend.text = element_text(size = 24),
      legend.title = element_text(size = 24)
    ),
  # list of plots
  labels = "auto",
  # labels
  common.legend = T,
  # common legend
  legend = "right",
  # legend position
  align = "hv",
  # Align them both, horizontal and vertical
  nrow = 2,
  # number of rows
  font.label = list(size = 24),
  ncol = 2  # number of columns
)

P7.allplot <- annotate_figure(P7allp,
                              top = text_grob("PORCE MODEL",
                                              color = "black",
                                              size = 28))
# Final Figure A2 - Distributions of age for each pool in Porce model and
# 14C distributions of pools in 1965, 2027 and 2100
P7.allplot

### Figure A3
# transit time distribution
TTP7dist <- ggplot(data = TTP7dens,
                   mapping = aes(x = age)) +
  geom_line(mapping = aes(y = TTP7.transitTimeDensity, 
                          color = "TT"),
            size = 1.3) +
  coord_cartesian(xlim = c(0, 150),
                  ylim = c(0, 0.1)) +
  labs(x = "Age (years)",
       y = "Probability density",
       color = element_blank()) +
  scale_color_manual(labels = c("TT" = "Outflux"),
                     values = c("TT" = "#F8766D")) +
  theme(legend.position = c(0.75, 0.65))
#14C distributions
allTTP7RDC <-
  ggplot(mapping = aes(x = mid, 
                       y = M)) +
  geom_bar(
    data = P7TThist1965,
    mapping = aes(color = "y1",
                  fill = "y1"),
    stat = "identity",
    position = position_dodge(),
    alpha = .4
  ) +
  geom_bar(
    data = P7TThist2027.p,
    mapping = aes(color = "y2",
                  fill = "y2"),
    stat = "identity",
    position = position_dodge(),
    alpha = .3
  ) +
  geom_bar(
    data = P7TThist2100.p,
    mapping = aes(color = "y3",
                  fill = "y3"),
    stat = "identity",
    position = position_dodge(),
    alpha = .3
  ) +
  labs(y = expression(paste("Mass of Carbon (Mg h",a^-1,")")),
       x = expression(paste(Delta ^ {14}, "C (\u2030)"))) +
  scale_color_manual(
    name = "Year of Observation",
    labels = c("y1" = "1965 AD",
               "y2" = "2027 AD",
               "y3" = "2100 AD"),
    values = c("y1" = "blue",
               "y2" = "darkolivegreen3",
               "y3" = "deeppink3")
  )  +
  scale_fill_manual(
    name = "Year of Observation",
    labels = c("y1" = "1965 AD",
               "y2" = "2027 AD",
               "y3" = "2100 AD"),
    values = c("y1" = "blue",
               "y2" = "darkolivegreen3",
               "y3" = "deeppink3")
  ) +
  theme(legend.position = c(0.8, 0.65))

# arranging alltogether
P7TTallp <- ggpubr::ggarrange(
  TTP7dist +
    ggtitle("Transit time distribution", subtitle = "Any year of observation") +
    theme(
      plot.title = element_text(hjust = 0.5,
                                size = 24,
                                face = "bold"),
      legend.position = "none",
      axis.title = element_text(size = 28),
      axis.text = element_text(size = 28),
      plot.subtitle = element_text(size = 24)
    ),
  allTTHFSRDC +
    ggtitle("Outflux radiocarbon distribution") +
    theme(
      plot.title = element_text(hjust = 0.5,
                                size = 24,
                                face = "bold"),
      legend.position = "top",
      axis.title = element_text(size = 28),
      axis.text = element_text(size = 28),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16)
    ),
  # list of plots
  labels = "auto",
  # labels
  common.legend = F,
  # common legend
  align = "v",
  # Align them both, horizontal and vertical
  nrow = 1,
  # number of rows
  font.label = list(size = 24),
  ncol = 2  # number of columns
)

P7TT.all <- annotate_figure(P7TTallp,
                            top = text_grob("PORCE MODEL",
                                            color = "black",
                                            size = 28))
# Final Figure A3 - Porce model transit time distribution and 14C distributions
# for the years 1965, 2027 and 2100
P7TT.all

### Figure A5
# age distributions pools
poolEmages <-
  ggplot(data = SAEmpool,
         mapping = aes(x = age)) +
  geom_line(mapping = aes(y = X1, 
                          color = "P1"),
            size = 1.1) +
  geom_line(mapping = aes(y = X2, 
                          color = "P2"),
            size = 1.1) +
  geom_line(mapping = aes(y = X3,
                          color = "P3"),
            size = 1.1) +
  geom_line(mapping = aes(y = X4, 
                          color = "P4"),
            size = 1.1) +
  geom_line(mapping = aes(y = X5,
                          color = "P5"),
            size = 1.1) +
  coord_cartesian(xlim = c(0, 150),
                  ylim = c(0, 0.3)) +
  labs(x = "Age (years)",
       y = "Probability density",
       color = "Pools") +
  scale_color_manual(
    labels = c(
      "P1" = bquote(paste(x[1], ": Non-woody Tree Parts")),
      "P2" = bquote(paste(x[2], ": Woody Tree Parts")),
      "P3" = bquote(paste(x[3], ": Ground Vegetation")),
      "P4" = bquote(paste(x[4], ": Detritus/Decomposers")),
      "P5" = bquote(paste(x[5], ": Active Soil Carbon"))
    ),
    values = c(
      "P1" = "darkolivegreen3",
      "P2" = "blue",
      "P3" = "#009E73",
      "P4" = "#F0E442",
      "P5" = "deeppink3"
    )
  ) +
  theme(legend.position = c(0.75, 0.65))
# 14C distributions pools
# AD 1965
PoolEm_RDC1965 <- ggplot(mapping = aes(x = mid,
                                       y = M)) +
  geom_bar(
    data = EmP1hist1965,
    mapping = aes(fill = "P1"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = EmP2hist1965,
    mapping = aes(fill = "P2"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = EmP3hist1965,
    mapping = aes(fill = "P3"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = EmP4hist1965,
    mapping = aes(fill = "P4"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = EmP5hist1965,
    mapping = aes(fill = "P5"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  )

LegPoolEm1965 <- PoolEm_RDC1965 +
  labs(y = "Mass of Carbon (Pg) ",
       x = expression(paste(Delta ^ {14}, "C (\u2030)")),
       fill = "Pools") +
  scale_fill_manual(
    labels = c(
      "P1" = poolnamesEm[1],
      "P2" = poolnamesEm[2],
      "P3" = poolnamesEm[3],
      "P4" = poolnamesEm[4],
      "P5" = poolnamesEm[5]
    ),
    values = c(
      "P1" = "darkolivegreen3",
      "P2" = "blue",
      "P3" = "#009E73",
      "P4" = "#F0E442",
      "P5" = "deeppink3"
    )
  ) +
  ggtitle("Year of observation = 1965") +
  theme(legend.position = c(0.75, 0.65))
# AD 2027
PoolEm_RDC2027 <- ggplot(mapping = aes(x = mid,
                                       y = M)) +
  geom_bar(
    data = EmP1hist2027,
    mapping = aes(fill = "P1"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = EmP2hist2027,
    mapping = aes(fill = "P2"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = EmP3hist2027,
    mapping = aes(fill = "P3"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = EmP4hist2027,
    mapping = aes(fill = "P4"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = EmP5hist2027,
    mapping = aes(fill = "P5"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  )

LegPoolEm2027 <- PoolEm_RDC2027 +
  labs(y = "Mass of Carbon (Pg) ",
       x = expression(paste(Delta ^ {14}, "C (\u2030)")),
       fill = "Pools") +
  scale_fill_manual(
    labels = c(
      "P1" = poolnamesEm[1],
      "P2" = poolnamesEm[2],
      "P3" = poolnamesEm[3],
      "P4" = poolnamesEm[4],
      "P5" = poolnamesEm[5]
    ),
    values = c(
      "P1" = "darkolivegreen3",
      "P2" = "blue",
      "P3" = "#009E73",
      "P4" = "#F0E442",
      "P5" = "deeppink3"
    )
  ) +
  ggtitle("Year of observation = 2027") +
  theme(legend.position = c(0.75, 0.65))
# AD 2100
PoolEm_RDC2100 <- ggplot(mapping = aes(x = mid,
                                       y = M)) +
  geom_bar(
    data = EmP1hist2100,
    mapping = aes(fill = "P1"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = EmP2hist2100,
    mapping = aes(fill = "P2"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = EmP3hist2100,
    mapping = aes(fill = "P3"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = EmP4hist2100,
    mapping = aes(fill = "P4"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  ) +
  geom_bar(
    data = EmP5hist2100,
    mapping = aes(fill = "P5"),
    stat = "identity",
    position = position_dodge(),
    color = "black"
  )

LegPoolEm2100 <- PoolEm_RDC2100 +
  labs(y = "Mass of Carbon (Pg) ",
       x = expression(paste(Delta ^ {14}, "C (\u2030)")),
       fill = "Pools") +
  scale_fill_manual(
    labels = c(
      "P1" = poolnamesEm[1],
      "P2" = poolnamesEm[2],
      "P3" = poolnamesEm[3],
      "P4" = poolnamesEm[4],
      "P5" = poolnamesEm[5]
    ),
    values = c(
      "P1" = "darkolivegreen3",
      "P2" = "blue",
      "P3" = "#009E73",
      "P4" = "#F0E442",
      "P5" = "deeppink3"
    )
  ) +
  ggtitle("Year of observation = 2100") +
  theme(legend.position = c(0.75, 0.65))
# arranging alltogether
Emallp <- ggpubr::ggarrange(
  poolEmages +
    ggtitle("Any year of observation") +
    theme(
      plot.title = element_text(hjust = 0.5,
                                size = 24,
                                face = "bold"),
      axis.title = element_text(size = 28),
      axis.text = element_text(size = 28),
      legend.text = element_text(size = 24),
      legend.title = element_text(size = 24)
    ),
  LegPoolEm1965 +
    ggtitle("Year of observation = 1965") +
    theme(
      plot.title = element_text(hjust = 0.5,
                                size = 24,
                                face = "bold"),
      axis.title = element_text(size = 28),
      axis.text = element_text(size = 28),
      legend.text = element_text(size = 24),
      legend.title = element_text(size = 24)
    ),
  LegPoolEm2027 +
    ggtitle("Year of observation = 2027") +
    theme(
      plot.title = element_text(hjust = 0.5,
                                size = 24,
                                face = "bold"),
      axis.title = element_text(size = 28),
      axis.text = element_text(size = 28),
      legend.text = element_text(size = 24),
      legend.title = element_text(size = 24)
    ),
  LegPoolEm2100 +
    ggtitle("Year of observation = 2100") +
    theme(
      plot.title = element_text(hjust = 0.5,
                                size = 24,
                                face = "bold"),
      axis.title = element_text(size = 28),
      axis.text = element_text(size = 28),
      legend.text = element_text(size = 24),
      legend.title = element_text(size = 24)
    ),
  # list of plots
  labels = "auto",
  # labels
  common.legend = T,
  # common legend
  legend = "right",
  # legend position
  align = "hv",
  # Align them both, horizontal and vertical
  nrow = 2,
  # number of rows
  ncol = 2,  # number of columns
  font.label = list(size = 24)) 

Em.allplot <- annotate_figure(Emallp,
                              top = text_grob(
                                "EMANUEL MODEL",
                                color = "black",
                                size = 28
                              ))
# Final Figure A5 - Distributions of age for each pool in Emanuel model and
# 14C distributions of pools in 1965, 2027 and 2100
Em.allplot

### Figure A6
# transit time distribution
TTEmdist <- ggplot(data = TTEmdens,
                   mapping = aes(x = age)) +
  geom_line(mapping = aes(y = TTEm.transitTimeDensity,
                          color = "TT"),
            size = 1.3) +
  coord_cartesian(xlim = c(0, 150),
                  ylim = c(0, 0.1)) +
  labs(x = "Age (years)",
       y = "Probability density",
       color = element_blank()) +
  scale_color_manual(
    labels = c("TT" = "Outflux"),
    values = c("TT" = "#F8766D")
  ) +
  theme(legend.position = c(0.75, 0.65))
#14C distributions
allTTEmRDC <-
  ggplot(mapping = aes(x = mid,
                       y = M)) +
  geom_bar(
    data = EmTThist1965,
    mapping = aes(color = "y1",
                  fill = "y1"),
    stat = "identity",
    position = position_dodge(),
    alpha = .4
  ) +
  geom_bar(
    data = EmTThist2027.p,
    mapping = aes(color = "y2",
                  fill = "y2"),
    stat = "identity",
    position = position_dodge(),
    alpha = .3
  ) +
  geom_bar(
    data = EmTThist2100.p,
    mapping = aes(color = "y3",
                  fill = "y3"),
    stat = "identity",
    position = position_dodge(),
    alpha = .3
  ) +
  labs(y = "Mass of Carbon (Pg)",
       x = expression(paste(Delta ^ {14}, "C (\u2030)"))) +
  scale_color_manual(
    name = "Year of Observation",
    labels = c("y1" = "1965 AD",
               "y2" = "2027 AD",
               "y3" = "2100 AD"),
    values = c("y1" = "blue",
               "y2" = "darkolivegreen3",
               "y3" = "deeppink3")
  )  +
  scale_fill_manual(
    name = "Year of Observation",
    labels = c("y1" = "1965 AD",
               "y2" = "2027 AD",
               "y3" = "2100 AD"),
    values = c("y1" = "blue",
               "y2" = "darkolivegreen3",
               "y3" = "deeppink3")
  ) +
  theme(legend.position = c(0.8, 0.65))
# arranging alltogether
EmTTallp <- ggpubr::ggarrange(
  TTEmdist +
    ggtitle("Transit time distribution", 
            subtitle = "Any year of observation") +
    theme(
      plot.title = element_text(hjust = 0.5,
                                size = 24,
                                face = "bold"),
      legend.position = "none",
      axis.title = element_text(size = 28),
      axis.text = element_text(size = 28),
      plot.subtitle = element_text(size = 24)
    ),
  allTTEmRDC +
    ggtitle("Outflux radiocarbon distribution") +
    theme(
      plot.title = element_text(hjust = 0.5,
                                size = 24,
                                face = "bold"),
      legend.position = "top",
      axis.title = element_text(size = 28),
      axis.text = element_text(size = 28),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16)
    ),
  # list of plots
  labels = "auto",
  # labels
  common.legend = F,
  # common legend
  align = "v",
  # Align them both, horizontal and vertical
  nrow = 1,
  # number of rows
  ncol = 2,
  # number of columns
  font.label = list(size = 24)
) 

EmTT.all <- annotate_figure(EmTTallp,
                            top = text_grob(
                              "EMANUEL MODEL",
                              color = "black",
                              size = 28))

# Final Figure A6 - Emanuel model transit time distribution and 14C distributions
# for the years 1965, 2027 and 2100
EmTT.all


