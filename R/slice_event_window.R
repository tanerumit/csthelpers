#' Slice Fixed-Length Event-Based Windows from a Timeline
#'
#' Constructs fixed-duration weather/climate windows anchored to a specified
#' start month and aligned with the year associated with each event sequence.
#' The window always spans exactly `total.days` days and is guaranteed to fit
#' fully within the provided `timeline`. This is useful for generating
#' scenario-based weather slices (e.g., for hydrology, climate stress-testing,
#' drought analysis, or impact models) where each slice begins on the same
#' calendar boundary (e.g., January 1, June 1) regardless of actual event dates.
#'
#' @description
#' For each event set:
#' \enumerate{
#'   \item Identify the last event date.
#'   \item Subtract `days.before` to estimate a reference year.
#'   \item Clamp this year to the range of valid years in which a full
#'         `total.days` window can fit inside `timeline`.
#'   \item Construct a window starting at `month.start`-01 of that year.
#'   \item Extend the window for `total.days` days.
#'   \item Return the indices of `timeline` that fall inside that window.
#' }
#' All windows are forced to equal length. If a single event vector is
#' provided, a single index vector is returned. If a list of event vectors
#' is provided, a list of index vectors is returned.
#'
#' @param events A `Date` vector or list of `Date` vectors. Each vector
#'   represents one event sequence (e.g., drought period, extreme spell).
#'
#' @param timeline A `Date` vector defining the complete historical range
#'   from which windows are extracted.
#'
#' @param days.before Integer. Number of days before the last event used to
#'   estimate the target year from which the slice will start.
#'   (This does *not* affect window length.)
#'
#' @param total.days Integer. Total length of the extracted window. The output
#'   window is always exactly this long, assuming `timeline` is long enough.
#'
#' @param month.start Integer (1–12). Month of the year on which each window
#'   starts (always day 1 of that month). For example, `1` = January 1,
#'   `6` = June 1.
#'
#' @return A vector of indices (if a single event is provided) or a list of
#'   equal-length index vectors (if multiple event sets are provided).
#'
#' @examples
#' # Example: extract 5-year windows starting January 1
#' ev  <- as.Date(c("1980-03-10", "1980-06-12"))
#' tl  <- seq.Date(as.Date("1980-01-01"), as.Date("2010-12-31"), by = "day")
#'
#' idx <- slice_event_window(
#'   events      = ev,
#'   timeline    = tl,
#'   days.before = 2 * 365,
#'   total.days  = 5 * 365,
#'   month.start = 1
#' )
#'
#' # Access the actual slice
#' tl[idx]
#'
#' @export
slice_event_window <- function(
    events,
    timeline,
    days.before = 365,
    total.days  = 5 * 365,
    month.start = 1
) {
  timeline <- as.Date(timeline)
  t_min <- min(timeline)
  t_max <- max(timeline)
  
  # Ensure timeline can support the required fixed-length window
  if (as.integer(t_max - t_min + 1) < total.days) {
    stop("timeline is shorter than total.days; cannot construct full-length windows.")
  }
  
  # Normalize to list
  if (!is.list(events)) {
    events <- list(events)
    single_input <- TRUE
  } else {
    single_input <- FALSE
  }
  events <- lapply(events, function(x) as.Date(sort(unique(x))))
  
  # Determine allowed start years
  earliest_year_allowed <- as.integer(format(t_min, "%Y"))
  earliest_start_date <- as.Date(paste0(earliest_year_allowed, "-", month.start, "-01"))
  if (earliest_start_date < t_min) {
    earliest_year_allowed <- earliest_year_allowed + 1
  }
  
  latest_year_allowed <- as.integer(format(t_max - (total.days - 1), "%Y"))
  latest_start_date <- as.Date(paste0(latest_year_allowed, "-", month.start, "-01"))
  if (latest_start_date > (t_max - (total.days - 1))) {
    latest_year_allowed <- latest_year_allowed - 1
  }
  
  if (earliest_year_allowed > latest_year_allowed) {
    stop("No valid start year exists that can accommodate total.days with the specified month.start.")
  }
  
  # Single-window generator
  get_window <- function(ev) {
    end_event <- max(ev)
    
    # Estimate target year
    estimate_start <- end_event - days.before
    target_year <- as.integer(format(estimate_start, "%Y"))
    
    # Clamp into valid range
    start_year <- max(min(target_year, latest_year_allowed), earliest_year_allowed)
    
    # Construct window start
    start_slice <- as.Date(sprintf("%04d-%02d-01", start_year, month.start))
    end_slice   <- start_slice + (total.days - 1)
    
    if (start_slice < t_min || end_slice > t_max) {
      stop("Internal error: constructed slice outside timeline.")
    }
    
    which(timeline >= start_slice & timeline <= end_slice)
  }
  
  raw <- lapply(events, get_window)
  
  # Equalize lengths
  max_len <- max(lengths(raw))
  out <- lapply(raw, function(idx) {
    if (length(idx) < max_len) {
      c(idx, rep(idx[length(idx)], max_len - length(idx)))
    } else {
      idx
    }
  })
  
  if (single_input) out <- out[[1]]
  out
}
