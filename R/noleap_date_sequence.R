#' @title Generate a No-Leap Daily Date Sequence
#'
#' @description
#' Creates a continuous sequence of dates that excludes February 29.
#' The sequence begins at a specified start date and continues for a
#' user-defined number of days. The function cycles through a
#' non-leap-year month-day template (Jan 1–Dec 31 without Feb 29) and
#' assigns appropriate year offsets.
#'
#' @param start_date Date. First day of the sequence.
#' @param n_days Integer. Number of days to generate after `start_date`.
#'
#' @return A `Date` vector of length `n_days` with no Feb 29 included.
#'
#' @details
#' The function constructs a month-day template from a non-leap year
#' (2001), removes Feb 29 for robustness, and then repeats this
#' template as needed to cover `n_days`. Year offsets are applied based
#' on the position in the template cycle.
#'
#' @examples
#' noleap_date_sequence("1980-01-01", 400)
#'
#' @export
noleap_date_sequence <- function(start_date, n_days) {
  # Ensure Date class
  start_date <- as.Date(start_date)
  
  # Build non-leap template year
  template <- seq.Date(as.Date("2001-01-01"), as.Date("2001-12-31"), by = "day")
  template <- template[format(template, "%m-%d") != "02-29"]
  
  # Extract month-day patterns
  month_day <- format(template, "%m-%d")
  
  # Repeat template to required length
  needed <- ceiling(n_days / length(month_day))
  md_seq <- rep(month_day, needed)[seq_len(n_days)]
  
  # Compute year offsets
  start_year <- as.integer(format(start_date, "%Y"))
  year_offsets <- as.integer((seq_len(n_days) - 1) / length(month_day))
  years <- start_year + year_offsets
  
  # Build final dates
  out <- as.Date(paste0(years, "-", md_seq))
  
  out
}
