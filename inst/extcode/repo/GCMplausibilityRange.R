#' Compute Plausibility Ranges for GCM Scenarios
#'
#' @description
#' Computes lower and upper bounds of a response surface for multiple GCM
#' scenarios and locations. The function constructs interpolated response
#' surfaces, extracts scenario-specific ellipses in \code{prcp–tavg} space,
#' identifies intersecting points, and evaluates the corresponding response
#' values to derive plausibility ranges.
#'
#' @param str.data Optional. Data frame containing historical or reference
#'   response-surface values used for interpolation. Structure must be
#'   compatible with the response variable extracted from \code{str_dfl}.
#' @param gcm.data Data frame of GCM outputs containing at least the variables
#'   \code{scenario}, \code{prcp}, and \code{tavg}. Used to compute scenario-
#'   specific ellipses.
#' @param gcm.scenarios Character vector of scenario identifiers to evaluate.
#' @param clevel Numeric. Confidence level for the ellipse constructed via
#'   \code{stat_ellipse()}. Default is \code{0.95}.
#' @param location.list Vector of location identifiers over which plausibility
#'   ranges are computed.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Generates a base grid of \code{prcp} and \code{tavg} values using
#'         interpolation over unique precipitation and temperature axes.
#'   \item Extracts points from this grid and determines which points lie inside
#'         the scenario-specific ellipses derived from \code{gcm.data}.
#'   \item For each location in \code{location.list}, constructs an interpolated
#'         response surface and evaluates the response variable at the
#'         ellipse-intersection points.
#'   \item Computes lower and upper bounds of the resulting response values.
#' }
#'
#' The function returns a data frame with one row per combination of scenario
#' and location.
#'
#' @return
#' A data frame containing:
#' \describe{
#'   \item{scenario}{Scenario identifier.}
#'   \item{location}{Location identifier.}
#'   \item{lo}{Lower plausibility bound.}
#'   \item{up}{Upper plausibility bound.}
#' }
#'
#' @note
#' This function relies on external objects such as \code{tavg_unique},
#' \code{prcp_unique}, and \code{str_dfl}. These must exist in the calling
#' environment. Ellipse extraction uses ggplot2 internals via
#' \code{ggplot_build()}.
#'
#' @examples
#' \dontrun{
#' GCMplausibilityRange(
#'   str.data = ref_surface,
#'   gcm.data = gcm_outputs,
#'   gcm.scenarios = c("ssp126", "ssp585"),
#'   clevel = 0.95,
#'   location.list = c("loc1", "loc2")
#' )
#' }
#'
#' @export
GCMplausibilityRange <- function(
    str.data = NULL,
    gcm.data = NULL,
    gcm.scenarios = NULL,
    clevel = 0.95,
    location.list ) {

  # Suppress warnings
  options(warn = -1)

  # Do for one location, metric and each gcm scenario
  OUTDF <- expand_grid(scenario = gcm.scenarios, location = location.list) %>%
    mutate(lo = NA, up = NA) %>%
    arrange(location)

  # EXTRACT CON VARIABLE ONCE ####################################################

  ### Do this for single location!
  con <- setNames(vector("list", length = length(gcm.scenarios)), gcm.scenarios)

  zmat <- matrix(0, nrow = length(tavg_unique), ncol = length(prcp_unique))
  fit  <- list(x = prcp_unique, y = tavg_unique, z = zmat)
  str0_i <- fields::interp.surface.grid(fit,
                                        grid.list = list(x = seq(min(prcp_unique), max(prcp_unique), length = 100),
                                                         y = seq(min(tavg_unique), max(tavg_unique), length = 100)))

  str0 <- tidyr::expand_grid(prcp = str0_i$x, tavg = str0_i$y)

  # Extract points from the response surface
  p <- ggplot(str0, aes(x = prcp, y = tavg)) +
    geom_point(aes())
  points <- ggplot_build(p)$data[[1]]

  # LOOP THROUGH GCM SCENARIOS
  for (s in 1:length(gcm.scenarios)) {

    gcm_data_s <- gcm.data %>% filter(scenario == gcm.scenarios[s])

    # Extract ellipse points from GCMs
    p1 <- p + stat_ellipse(data = gcm_data_s, aes(x = prcp, y = tavg), level = clevel, type = "norm", color = "red", size = 2)
    ell <- ggplot_build(p1)$data[[2]]

    # Find intersecting points on the ellipse
    con[[gcm.scenarios[s]]] <- which(as.logical(sp::point.in.polygon(points$x, points$y, ell$x, ell$y)))

  }

  ######### LOOP THROUGH EACH LOCATION
  for (x in 1:nrow(OUTDF)) {

    # Current identifiers
    sx <- OUTDF$scenario[x]
    lx <- OUTDF$location[x]

    # Interpolate results for the current location
    z <- str_dfl %>% filter(location == lx) %>% pull(value)

    zmat <- matrix(z, nrow = length(tavg_unique),
                   ncol = length(prcp_unique), byrow = TRUE)
    fit <- list(x = prcp_unique, y = tavg_unique, z = zmat)

    str_dfl_interpi <- fields::interp.surface.grid(fit,
                                                   grid.list = list(x = seq(min(prcp_unique), max(prcp_unique), length = 100),
                                                                    y = seq(min(tavg_unique), max(tavg_unique), length = 100)))

    str_dfl_interp <- tidyr::expand_grid(prcp = str_dfl_interpi$x, tavg = str_dfl_interpi$y) %>%
      mutate(z = as.vector(str_dfl_interpi$z))

    # Find matching values
    valrange <- str_dfl_interp[con[[sx]], ]$z
    OUTDF$lo[x] <- min(valrange)
    OUTDF$up[x] <- max(valrange)

  }

  OUTDF

}

