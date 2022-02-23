#' TTRDC
#'
#' @param A compartmental matrix
#' @param u input vector
#' @param nyears number of years
#' @param yrobs year of observation
#' @param C14atm atmospheric radiocarbon content
#' @param C14YearAD calendar year (in AD) corresponding to C14atm
#' @param h discretization size
#'
#' @return
#' @export

if (!require("SoilR")) install.packages("SoilR")
if (!require("FME")) install.packages("FME")

library(SoilR)
library(FME)

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
  D14C.decay <- ((((newResC14curve / 1000) + 1) * exp((-1 / 8267) * tt)) - 1) * 1000 
  return(data.frame(YearAD = year.ad,
                    D14C = D14C.decay,
                    Mass = TTAMD * h))
}