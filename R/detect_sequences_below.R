

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
