
#' Calculate the Standardized Precipitation Evapotranspiration Index (SPEI)
#'
#' @description
#' Computes the Standardized Precipitation Evapotranspiration Index (SPEI)
#' from monthly precipitation and potential evapotranspiration (PET) data,
#' using the statistical framework implemented in the \pkg{SPEI} package.
#' 
#' The function first computes the climatic water balance \eqn{D = P - PET},
#' converts it into a monthly time series, and applies a probability
#' distribution (by default, log-logistic) to obtain standardized
#' anomalies of moisture deficit or surplus at a chosen accumulation scale.
#'
#' @param P Numeric vector of monthly precipitation values (mm).
#' @param PET Numeric vector of monthly potential evapotranspiration values (mm).
#' @param scale Integer. Accumulation period in months (e.g., 1, 3, 6, 12).
#'   This defines the temporal scale over which the water balance is aggregated
#'   before standardization.
#' @param dist Character string indicating the distribution to fit.
#'   Common choices are `"log-Logistic"` (default) or `"Gamma"`.
#'   See \link[SPEI:spei]{SPEI::spei()} for available options.
#' @param ref_start Optional numeric vector of length two (e.g., \code{c(1981, 1)})
#'   specifying the start year and month of the reference period used for
#'   fitting the distribution. If \code{NULL}, the entire record is used.
#' @param ref_end Optional numeric vector of length two (e.g., \code{c(2010, 12)})
#'   specifying the end year and month of the reference period used for
#'   fitting the distribution.
#' @param dates Optional vector of class \code{Date} corresponding to each
#'   observation. If provided, the function automatically determines the
#'   start year and month for the internal time-series object. If omitted,
#'   the series is assumed to begin at (1, 1) with monthly frequency.
#'
#' @details
#' The function acts as a simplified wrapper around \code{\link[SPEI:spei]{SPEI::spei()}}.
#' It is designed for convenience in workflows where precipitation (P),
#' evapotranspiration (PET), and dates are stored as plain vectors rather
#' than time-series objects. The resulting SPEI values are standardized
#' z-scores representing the relative moisture deficit or surplus at each time step.
#'
#' Positive SPEI values indicate wetter-than-normal conditions,
#' while negative values indicate drier-than-normal (potential drought)
#' conditions.
#'
#' @return
#' A numeric vector of SPEI values (same length as input),
#' corresponding to the fitted index at each time step.
#'
#' @references
#' Beguería, S., Vicente-Serrano, S. M., Reig, F., & Latorre, B. (2014).
#' Standardized precipitation evapotranspiration index (SPEI) revisited:
#' parameter fitting, evapotranspiration models, tools, datasets and trends worldwide.
#' \emph{International Journal of Climatology, 34}(10), 3001–3023.
#'
#' Vicente-Serrano, S. M., Beguería, S., & López-Moreno, J. I. (2010).
#' A multiscalar drought index sensitive to global warming:
#' The Standardized Precipitation Evapotranspiration Index.
#' \emph{Journal of Climate, 23}(7), 1696–1718.
#'
#' @seealso
#' \code{\link[SPEI:spei]{SPEI::spei()}} for the underlying implementation
#' and \code{\link[SPEI:spi]{SPEI::spi()}} for the precipitation-only variant.
#'
#' @examples
#' \dontrun{
#' library(SPEI)
#' set.seed(42)
#' dates <- seq(as.Date("1981-01-01"), as.Date("2020-12-01"), by = "month")
#' P   <- rnorm(length(dates), 60, 15)  # precipitation (mm)
#' PET <- rnorm(length(dates), 45, 10)  # PET (mm)
#'
#' spei6 <- calc_spei(P, PET, scale = 6, dates = dates,
#'                    ref_start = c(1981, 1), ref_end = c(2010, 12))
#'
#' plot(dates, spei6, type = "l", col = "steelblue",
#'      ylab = "SPEI (6-month)", xlab = "Year")
#' abline(h = c(-2, -1, 0, 1, 2), lty = 2, col = "gray")
#' }
#'
#' @export
calc_spei <- function(P, PET, scale = 3,
                      dist = "log-Logistic",
                      ref_start = NULL,
                      ref_end = NULL,
                      dates = NULL) {
  if (!requireNamespace("SPEI", quietly = TRUE))
    stop("Package 'SPEI' is required. Install with install.packages('SPEI').")
  
  stopifnot(length(P) == length(PET))
  D <- P - PET
  
  if (!is.null(dates)) {
    start_year  <- as.numeric(format(min(dates), "%Y"))
    start_month <- as.numeric(format(min(dates), "%m"))
    D <- ts(D, start = c(start_year, start_month), frequency = 12)
  } else {
    D <- ts(D, frequency = 12)
  }
  
  # Run SPEI silently (suppress cat, message, warning output)
  spei_obj <- suppressWarnings(
    suppressMessages(
      {
        tmp <- capture.output(
          res <- SPEI::spei(D,
                            scale = scale,
                            distribution = dist,
                            ref.start = ref_start,
                            ref.end = ref_end),
          file = NULL
        )
        res
      }
    )
  )
  
  out <- as.numeric(spei_obj$fitted)
  return(out)
}



