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
