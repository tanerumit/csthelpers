#' Find Nearest Neighbor (k = 1) in Euclidean Space
#'
#' @description
#' A simple, pure-R implementation of a single nearest-neighbor search
#' between two numeric coordinate sets.
#' For each query point, the function identifies the index (and optionally
#' the distance) of the closest point in the reference dataset.
#'
#' @param data Numeric matrix or data frame containing reference points,
#'   with one row per point and one column per coordinate dimension (e.g. x, y, z).
#' @param query Numeric matrix or data frame containing query points to match,
#'   with the same number of columns as `data`.
#' @param return_dist Logical; if `TRUE`, also returns the Euclidean distance
#'   to the nearest neighbor for each query point. Default is `FALSE`.
#'
#' @return
#' If `return_dist = FALSE`, returns an integer vector of indices indicating the
#' nearest neighbor in `data` for each query point.
#' If `return_dist = TRUE`, returns a list with two elements:
#' \describe{
#'   \item{idx}{Integer vector of nearest-neighbor indices.}
#'   \item{dist}{Numeric vector of Euclidean distances corresponding to each match.}
#' }
#'
#' @details
#' This function performs a brute-force nearest-neighbor search by computing
#' the squared Euclidean distance between each query point and all reference
#' points.
#' It is straightforward and accurate, but scales as O(n_query × n_data), so it
#' may be slow for very large datasets. For large-scale applications, consider
#' using [RANN::nn2()] or [FNN::get.knnx()] for optimized k-d tree implementations.
#'
#' @examples
#' set.seed(123)
#' data  <- data.frame(x = runif(10), y = runif(10))
#' query <- data.frame(x = c(0.2, 0.8), y = c(0.3, 0.9))
#'
#' # Return only indices of nearest neighbors
#' idx <- find_nn1(data, query)
#' idx
#'
#' # Return both indices and distances
#' res <- find_nn1(data, query, return_dist = TRUE)
#' res
#'
#' @seealso
#' [RANN::nn2()] and [FNN::get.knnx()] for faster approximate or exact nearest-neighbor search.
#'
#' @export
find_nn1 <- function(data, query, return_dist = FALSE) {
  data_mat  <- as.matrix(data)
  query_mat <- as.matrix(query)

  if (ncol(data_mat) != ncol(query_mat))
    stop("data and query must have the same number of columns.")

  n_query <- nrow(query_mat)
  n_data  <- nrow(data_mat)

  idx  <- integer(n_query)
  dist <- numeric(n_query)

  # brute-force nearest neighbor
  for (i in seq_len(n_query)) {
    # compute squared Euclidean distances
    d2 <- colSums((t(data_mat) - query_mat[i, ])^2)
    j <- which.min(d2)
    idx[i]  <- j
    dist[i] <- sqrt(d2[j])
  }

  if (return_dist)
    return(list(idx = idx, dist = dist))
  else
    return(idx)
}




#' Select a Space-Filling Subset of Points Using Maximin Sampling
#'
#' Implements a greedy maximin (farthest-point) selection algorithm to extract
#' a subset of \code{k} points from a two-dimensional point cloud. The algorithm
#' selects points such that each newly added point maximizes the minimum
#' Euclidean distance to all previously selected points, producing a
#' well-distributed, space-filling subset across the \code{(x, y)} domain.
#'
#' @description
#' Given two numeric vectors \code{x} and \code{y} of equal length, this function:
#' \enumerate{
#'   \item Forms a point set \eqn{(x_i, y_i)}.
#'   \item Computes the full pairwise distance matrix.
#'   \item Initializes the selection with the point closest to the centroid.
#'   \item Iteratively adds the point with the largest minimum distance from the
#'         already-selected set.
#' }
#' The result is a set of indices corresponding to the selected points. This
#' method is commonly used in experimental design, scenario reduction, and
#' space-filling sampling for climate and hydrological stress-testing.
#'
#' @param x Numeric vector. First coordinate (e.g., CMDmin).
#' @param y Numeric vector. Second coordinate (e.g., DPI).
#' @param k Integer. Number of points to select. Must satisfy \code{k <= length(x)}.
#'
#' @return
#' An integer vector of indices specifying the selected subset of points.
#'
#' @examples
#' set.seed(1)
#' x <- runif(50, 0, 100)
#' y <- runif(50, 0, 1)
#' idx <- select_maximin(x, y, k = 10)
#' plot(x, y, pch = 19)
#' points(x[idx], y[idx], col = "red", pch = 19, cex = 1.4)
#'
#' @export
select_maximin <- function(x, y, k) {
  pts <- cbind(x, y)
  n   <- nrow(pts)
  if (k > n) stop("k cannot be larger than number of points")

  # Precompute full distance matrix
  D <- as.matrix(dist(pts))

  # Start with a point near the centroid
  start <- which.min((x - mean(x))^2 + (y - mean(y))^2)
  selected  <- start
  remaining <- setdiff(seq_len(n), selected)

  while (length(selected) < k) {
    dmin <- apply(D[remaining, selected, drop = FALSE], 1, min)
    next_pt <- remaining[which.max(dmin)]
    selected  <- c(selected, next_pt)
    remaining <- setdiff(remaining, next_pt)
  }

  selected
}

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
#' Find Consecutive Sequences Below a Threshold
#'
#' Identifies all sequences of consecutive values in a numeric vector
#' that are below a specified threshold and have a minimum desired length.
#' Optionally, only sequences starting after a given index are considered.
#'
#' @param x Numeric vector to analyze.
#' @param threshold Numeric value. Elements strictly below this value are considered part of a sequence.
#' @param min_length Integer. Minimum number of consecutive elements below the threshold required for a sequence to be returned.
#' @param min_start_index Integer. Optional index threshold. Sequences starting before this index are ignored.
#' Default is \code{1} (i.e., all sequences considered).
#'
#' @return A list of integer vectors. Each element corresponds to the indices of one qualifying sequence in the original vector.
#' If no qualifying sequences are found, returns an empty list.
#'
#' @examples
#' x <- c(rep(5, 1000), 2, 1, 1, 0, 5, 2, 2, 2, 2)
#' find_sequences_below(x, threshold = 3, min_length = 3, min_start_index = 1000)
#'
#' @export
find_sequences_below <- function(x, threshold, min_length = 3, min_start_index = 1) {
  # Convert to logical vector: TRUE where x < threshold
  below <- x < threshold

  # Get run-length encoding of TRUE/FALSE segments
  rle_x <- rle(below)

  # Identify runs that are TRUE and long enough
  long_runs <- which(rle_x$values & rle_x$lengths >= min_length)

  # Convert RLE positions back to original indices
  if (length(long_runs) == 0) return(list())

  ends <- cumsum(rle_x$lengths)
  starts <- ends - rle_x$lengths + 1

  # Build list of index sequences
  result <- lapply(long_runs, function(i) seq(starts[i], ends[i]))

  # Filter by min_start_index safely
  keep <- vapply(result, function(idx) min(idx) >= min_start_index, logical(1))
  result <- result[keep]

  # Return empty list if nothing passes filter
  if (length(result) == 0) return(list())

  return(result)
}
