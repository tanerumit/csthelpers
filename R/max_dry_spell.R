#' Compute Maximum Dry-Spell Duration
#'
#' Calculates the longest continuous sequence of dry days in a precipitation
#' time series. A dry day is defined as any day with precipitation below a user-
#' specified threshold. The function uses run-length encoding to identify
#' consecutive dry runs and returns the maximum duration in days.
#'
#' @param P Numeric vector of daily precipitation values.
#' @param threshold Numeric. Precipitation threshold (in mm) used to define a
#'   dry day. Defaults to `1` (i.e., days with P < 1 mm are considered dry).
#'
#' @return
#' An integer value representing the length (in days) of the longest continuous
#' dry spell. Returns `0` if no dry days occur in the series.
#'
#' @examples
#' # Example precipitation series
#' P <- c(0, 0.2, 0.5, 3, 0, 0, 0, 5)
#'
#' # Longest dry spell with default threshold (1 mm)
#' max_dry_spell(P)
#'
#' # Using a stricter threshold
#' max_dry_spell(P, threshold = 0.1)
#'
#' @export
max_dry_spell <- function(P, threshold = 1) {

  dry <- P < threshold

  r <- rle(dry)

  dry_lengths <- r$lengths[r$values == TRUE]

  if (length(dry_lengths) == 0) return(0)

  max(dry_lengths)
  #quantile(dry_lengths, probs = 0.95)
}



compute_dpi <- function(spei_values, threshold = -1) {
  # Remove NAs
  x <- spei_values[!is.na(spei_values)]
  if (length(x) == 0) return(0)
  mean(x < threshold)
}



