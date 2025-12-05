#' Slice Time Windows Around Events
#'
#' Extracts backward-looking time windows around one or more sets of event dates,
#' optionally expands the window by a fixed number of days, and equalizes the
#' resulting index vectors to identical lengths. Each window is clipped to the
#' range of `full_dates`.
#'
#' @param event_dates A `Date` vector or a list of `Date` vectors. Each element
#'   represents a set of event dates. If a single vector is supplied, it is
#'   wrapped into a list internally.
#' @param full_dates A `Date` vector defining the complete timeline from which
#'   indices are selected.
#' @param window_years Numeric. Number of years to look back from the last event
#'   date in each set.
#' @param expand_days Integer. Optional symmetric padding (in days) added to
#'   both the start and end dates of each window after the year-based offset.
#'
#' @return If `event_dates` is a vector, returns an integer vector of indices
#'   corresponding to positions in `full_dates` that fall within the constructed
#'   window. If `event_dates` is a list, returns a list of such integer vectors,
#'   each padded to the same length.
#'
#' @details
#' For each set of event dates, the function:
#' \enumerate{
#'   \item Converts inputs to `Date`.
#'   \item Identifies the last event date.
#'   \item Creates a backward-looking window of `window_years` using a
#'     year-offset operation.
#'   \item Applies optional `expand_days` padding.
#'   \item Clips the window to the bounds of `full_dates`.
#'   \item Returns the indices of all `full_dates` inside the window.
#' }
#' When multiple event-date sets are provided, the index vectors are padded with
#' repetitions of their final index to ensure equal length.
#'
#' @importFrom lubridate `%m-%` years
#'
#' @examples
#' # Example event dates
#' ev  <- as.Date(c("2005-03-10", "2005-03-12"))
#' all <- seq.Date(as.Date("2000-01-01"), as.Date("2010-12-31"), by = "day")
#'
#' slice_time_window(ev, all, window_years = 3, expand_days = 10)
#'
#' @export
slice_time_window <- function(event_dates, full_dates,
                              window_years = 5, expand_days = 0) {
  # Ensure Date class
  full_dates  <- as.Date(full_dates)

  # Normalize: if event_dates is not a list, make it a list
  if (!is.list(event_dates)) {
    event_dates <- list(event_dates)
    single_input <- TRUE
  } else {
    single_input <- FALSE
  }

  # Convert all event date vectors to Date
  event_dates <- lapply(event_dates, as.Date)

  # ---- Core window extractor ----
  get_window <- function(ev_dates) {

    # End of the event block = last event date
    end_date <- max(ev_dates)

    # Backward-looking window
    start_date <- lubridate::`%m-%`(end_date, lubridate::years(window_years))

    # Apply optional padding
    start_date <- start_date - expand_days
    end_date   <- end_date   + expand_days

    # Clip to full_dates range
    start_date <- max(start_date, min(full_dates))
    end_date   <- min(end_date,   max(full_dates))

    # Return indices of full_dates in the window
    which(full_dates >= start_date & full_dates <= end_date)
  }

  # Extract raw windows
  raw_windows <- lapply(event_dates, get_window)

  # ---- Make all windows equal length ----
  max_len <- max(lengths(raw_windows))

  equalized <- lapply(raw_windows, function(idx) {
    n <- length(idx)
    if (n < max_len) {
      # pad with the *last index* to extend length
      c(idx, rep(idx[length(idx)], max_len - n))
    } else {
      idx
    }
  })

  if (single_input) equalized <- equalized[[1]]

  return(equalized)
}

