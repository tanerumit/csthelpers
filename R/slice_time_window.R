

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

