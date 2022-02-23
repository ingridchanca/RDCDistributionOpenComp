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

## Libraries
library(tidyverse)
library(SoilR) #for calculations of distributions of ages, transit time and 14C
library(FME) #for the calculation of distrib (works with SoilR)
library(Hmisc) #statistics (use for variance of theoretical 14C)
library(ggplot2) #for plots
library(grid) #arranging plots
library(gridExtra) #arranging plots
library(cowplot) #arranging plots

## directory ????
base <- file.path("~/RDCDistributionsOpenComp/")
data <- file.path(base, "data")

## plot theme used in the manuscript
theme_set(theme_classic(base_size = 18))

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

### Atmospheric Delta14CO2

#### Tropics -> IntCal20 + Graven et al. (2017), GMD = records for the tropics
#### + forecast RCP8.5 Graven (2015), PNAS
C14_Trpcs <- read.table(file = file.path(data, "C14Graven_Tropics_forecast_IntCal20")) # use for Porce model

#### Northern Hemisphere -> IntCal20 + Graven et al. (2017), GMD = records for NH
#### + forecast RCP8.5 Graven (2015), PNAS
C14_NH <- read.table(file = file.path(data, "C14Graven_PNAS_GMD_Int20_2100_-53050")) # use for HFS and Emanuel models





## Plots
### Figure 3 - Scheme of 14C curves for the period -53050 to 2100 and data sets used
nhtrpcs.curves <- data.frame(C14_Trpcs, C14_NH) # none=Tropics, 1=NH
C14nhtrpcs.p <- ggplot(data = nhtrpcs.curves) +
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

### Figure 4 - 

