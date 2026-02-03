
#------------------------------------------------------------------------------
# Shared helpers
#------------------------------------------------------------------------------

#' Median grid step for a (quasi-)regular axis
#' @keywords internal
# -------------------------------------------------------------------------
# Shared internal helpers: 2D KDE bandwidth selection and validation
# -------------------------------------------------------------------------


.kde_is_pd_2x2 <- function(M, tol = .Machine$double.eps) {
  if (!is.matrix(M) || any(dim(M) != c(2L, 2L))) return(FALSE)
  if (any(!is.finite(M))) return(FALSE)
  ev <- tryCatch(eigen(M, symmetric = TRUE, only.values = TRUE)$values,
                 error = function(e) NULL)
  if (is.null(ev) || any(!is.finite(ev))) return(FALSE)
  all(ev > tol)
}

.kde_try_Hpi_2d <- function(data_mat) {
  if (!requireNamespace("ks", quietly = TRUE)) return(NULL)
  if (!is.matrix(data_mat) || ncol(data_mat) != 2L) return(NULL)
  tryCatch({
    H <- ks::Hpi(x = data_mat, nstage = 2, pilot = "samse", pre = "sphere")
    if (!.kde_is_pd_2x2(H)) return(NULL)
    H
  }, error = function(e) NULL)
}

.kde_try_nrd_2d <- function(x, y) {
  if (!requireNamespace("MASS", quietly = TRUE)) return(NULL)
  tryCatch({
    bwx <- MASS::bandwidth.nrd(x)
    bwy <- MASS::bandwidth.nrd(y)
    if (!is.finite(bwx) || !is.finite(bwy) || bwx <= 0 || bwy <= 0) return(NULL)
    diag(c(bwx^2, bwy^2))
  }, error = function(e) NULL)
}

.kde_pick_H_2d <- function(data_mat,
                           bw_method = c("auto", "plugin", "nrd"),
                           bw_adjust = 1.0) {
  bw_method <- match.arg(bw_method)
  if (!is.finite(bw_adjust) || bw_adjust <= 0) return(NULL)

  H <- NULL
  if (bw_method %in% c("auto", "plugin")) {
    H <- .kde_try_Hpi_2d(data_mat)
  }
  if (is.null(H)) {
    H <- .kde_try_nrd_2d(data_mat[, 1], data_mat[, 2])
  }
  if (is.null(H)) return(NULL)

  H_adj <- H * (bw_adjust^2)
  if (!.kde_is_pd_2x2(H_adj)) return(NULL)
  H_adj
}
