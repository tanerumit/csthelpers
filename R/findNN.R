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
#' idx <- findNN(data, query)
#' idx
#'
#' # Return both indices and distances
#' res <- findNN(data, query, return_dist = TRUE)
#' res
#'
#' @seealso
#' [RANN::nn2()] and [FNN::get.knnx()] for faster approximate or exact nearest-neighbor search.
#'
#' @export
findNN <- function(data, query, return_dist = FALSE) {
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