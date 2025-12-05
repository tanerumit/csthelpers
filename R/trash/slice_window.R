#' Slice Daily Time Windows Before and After Event Periods
#'
#' Constructs fixed-length time windows around one or more event periods,
#' using a backward-looking window defined in days and an additional
#' forward-looking tail, also in days. The window is always anchored on the
#' *last* event date. All returned index vectors are padded to identical length
#' to ensure consistent downstream processing.
#'
#' @description
#' For each event-date set, the function:
#' \enumerate{
#'   \item Converts inputs to `Date`.
#'   \item Identifies the final (maximum) event date.
#'   \item Computes `start = end_event - days_before`.
#'   \item Computes `end   = end_event + days_after`.
#'   \item Clips both boundaries to the range of `timeline`.
#'   \item Returns the indices of `timeline` that lie within the window.
#' }
#' If multiple event-date sets are provided, output vectors are padded by
#' repeating the last index until all windows have equal length.
#'
#' @param events A `Date` vector or a list of `Date` vectors. Each element
#'   represents an event period. A single vector is automatically wrapped into
#'   a list.
#'
#' @param timeline A `Date` vector defining the full available daily sequence
#'   from which indices are selected.
#'
#' @param days_before Integer. Number of days to include *before* the last
#'   event date. Defines the backward-looking length of the window.
#'
#' @param days_after Integer. Number of days to include *after* the last
#'   event date. Adds a forward tail to the window.
#'
#' @return
#' If `events` is a single vector, returns one integer vector of indices into
#' `timeline`.
#' If `events` is a list, returns a list of such index vectors, all padded to
#' equal length.
#'
#' @examples
#' ev  <- as.Date(c("2002-01-01", "2002-02-01"))
#' tl  <- seq.Date(as.Date("1995-01-01"), as.Date("2005-12-31"), by = "day")
#'
#' # 3-year backward window (≈ 1095 days) + 1 year forward padding
#' out <- slice_window(ev, tl, days_before = 1095, days_after = 365)
#' length(out)
#'
#' @export
slice_window <- function(
    events,
    timeline,
    days_before = 365,
    days_after  = 0
) {
  timeline <- as.Date(timeline)

  # Normalize input to list
  if (!is.list(events)) {
    events <- list(events)
    single_input <- TRUE
  } else {
    single_input <- FALSE
  }

  # Standardize event-date sets
  events <- lapply(events, function(x) as.Date(sort(unique(x))))

  # Helper to clip dates
  clip_date <- function(x) {
    if (x < min(timeline)) return(min(timeline))
    if (x > max(timeline)) return(max(timeline))
    x
  }

  # Construct window for one event set
  get_window <- function(ev) {
    end_event   <- max(ev)
    start_slice <- end_event - days_before
    end_slice   <- end_event + days_after

    start_slice <- clip_date(start_slice)
    end_slice   <- clip_date(end_slice)

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
