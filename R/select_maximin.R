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

