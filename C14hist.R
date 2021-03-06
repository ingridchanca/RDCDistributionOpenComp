#' C14hist
#'
#' @param D14C Delta14C values obtained through PoolRDC, SystemRDC or TTRDC
#' @param Mass Corresponding mass densities from PoolRDC, SystemRDC or TTRDC
#' @param bin Bin size (b) of the histograms to be generated by the funtion
#'
#' @return
#' @export

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