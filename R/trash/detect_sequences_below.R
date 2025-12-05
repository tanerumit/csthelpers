#' Detect Sequences of Values Below a Threshold
#'
#' @description
#' Parses a character vector containing numeric values (possibly with R-style
#' index labels such as \code{"[17]"}), converts it to a numeric sequence, and
#' identifies all contiguous index sequences where values fall below a specified
#' threshold for at least a minimum duration.
#'
#' @param x Character vector. May include R index labels (e.g., \code{"[1]"}),
#'   whitespace, or multi-line printed numeric output. All values are extracted,
#'   cleaned, and converted to numeric.
#' @param threshold Numeric. Values strictly less than this threshold are
#'   considered part of a "below-threshold" run.
#' @param duration Integer. Minimum length of a contiguous below-threshold run
#'   required for inclusion in the output.
#'
#' @details
#' The function removes printed index labels, splits the cleaned character
#' vector into individual tokens, converts them to numeric values, and applies
#' run-length encoding to detect contiguous segments where values fall below the
#' threshold. Only segments with length greater than or equal to \code{duration}
#' are returned.
#'
#' @return
#' A list where each element is an integer vector of indices corresponding to a
#' contiguous run of values meeting the below-threshold and duration criteria.
#' If no runs qualify, an empty list is returned.
#'
#' @examples
#' x <- c("[1] 5 4 3 6 2 1 1 7")
#' detect_sequences_below(x, threshold = 3, duration = 2)
#'
#' @export
detect_sequences_below <- function(x, threshold, duration) {
  # Remove R index labels like "[17]" and collapse components
  x_clean <- gsub("\\[[0-9]+\\]", "", x)
  x_clean <- unlist(strsplit(x_clean, "\\s+"))
  x_clean <- x_clean[x_clean != ""]

  # Convert to numeric vector
  vals <- as.numeric(x_clean)

  # Identify positions below threshold
  below <- vals < threshold

  # Run-length encoding to detect sequences
  r <- rle(below)

  # Compute starting positions of each run
  run_starts <- cumsum(c(1, head(r$lengths, -1)))

  # Identify runs that meet both criteria
  qualifying <- which(r$values == TRUE & r$lengths >= duration)

  # Build list of index sequences
  result <- lapply(qualifying, function(i) {
    start <- run_starts[i]
    end   <- start + r$lengths[i] - 1
    seq(start, end)
  })

  return(result)
}
