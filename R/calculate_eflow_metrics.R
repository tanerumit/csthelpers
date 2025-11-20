#' Calculate Environmental Flow (E-flow) Hydrologic Metrics
#'
#' @description
#' Computes 70 eco-hydrological indicators describing magnitude, timing,
#' frequency, rate of change, and intermittency of daily discharge time-series.
#' Supports both single-site (vector) input and multi-site (matrix/data.frame)
#' input. Metrics are summarised per year and aggregated to median and IQR.
#'
#' @param Q Numeric vector, matrix, or data.frame. Daily discharge values.
#'   If matrix/data.frame: rows = dates, columns = sites.
#' @param date Vector of dates (`Date` or coercible to `Date`). Length must
#'   match `nrow(Q)` for multi-site input, or `length(Q)` for single-site.
#' @param high_thresh_quantile Upper quantile for high-flow threshold (default 0.9).
#' @param low_thresh_quantile Lower quantile for low-flow threshold (default 0.1).
#' @param zero_thr Lower quantile for zero-flow classification
#'   (default 0.02).
#'
#' @return
#' A tibble containing:
#' \itemize{
#'   \item \code{site} — site identifier (multi-site only)
#'   \item \code{metric_id} — numeric metric index (1–70)
#'   \item \code{group} — metric group (HF, MF, LF, RFC, IF)
#'   \item \code{metric} — short metric name
#'   \item \code{value} — computed value
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' dates <- seq.Date(as.Date("1990-01-01"), as.Date("2019-12-31"), by="day")
#' Q <- pmax(0, 40 + 10*sin(2*pi*(1:length(dates))/365) + rnorm(length(dates)))
#'
#' out <- calculate_eflow_metrics(Q, dates)
#' head(out)
#' }
#'
#' @import dplyr lubridate zoo tibble
#' @export
calculate_eflow_metrics <- function(
    Q, date,
    high_thresh_quantile = 0.9,
    low_thresh_quantile  = 0.1,
    zero_thr = 0.02
) {
  
  # --------------------------------------------------------------------------
  # INTERNAL FUNCTION: compute metrics for a single site
  # --------------------------------------------------------------------------
  calc_site <- function(Qv, datev) {
    
    stopifnot(length(Qv) == length(datev))
    datev <- as.Date(datev)
    
    # Precompute vectors
    Qv         <- as.numeric(Qv)
    yr         <- year(datev)
    mon        <- month(datev)
    jday       <- yday(datev)
    years      <- sort(unique(yr))
    ny         <- length(years)
    
    # Quantile thresholds
    high_thr   <- quantile(Qv, probs = high_thresh_quantile, na.rm = TRUE)
    low_thr    <- quantile(Qv, probs = low_thresh_quantile,  na.rm = TRUE)

    # Precompute rolling means once
    r3   <- zoo::rollmean(Qv,  3, align="right", fill=NA)
    r7   <- zoo::rollmean(Qv,  7, align="right", fill=NA)
    r30  <- zoo::rollmean(Qv, 30, align="right", fill=NA)
    r90  <- zoo::rollmean(Qv, 90, align="right", fill=NA)
    
    # Rate of change (same for all years)
    dQ <- diff(Qv)
    
    # Storage container
    res <- vector("list", ny)
    
    seq_stats <- function(cond) {
      r <- rle(cond)
      if (!any(r$values)) return(list(num=0, avg=0, start=NA))
      num <- sum(r$values)
      avg <- mean(r$lengths[r$values])
      # longest sequence start
      idx  <- which(r$values)[ which.max(r$lengths[r$values]) ]
      start <- sum(r$lengths[1:(idx-1)]) + 1
      list(num=num, avg=avg, start=start)
    }
    
    # ------------------------------------------------------------------------
    # YEAR LOOP
    # ------------------------------------------------------------------------
    for (i in seq_len(ny)) {
      
      y <- years[i]
      idx <- which(yr == y)
      if (length(idx) == 0) next
      
      Qy   <- Qv[idx]
      r3y  <- r3[idx]
      r7y  <- r7[idx]
      r30y <- r30[idx]
      r90y <- r90[idx]
      jdy  <- jday[idx]
      mony <- mon[idx]
      
      # high/low pulses
      hp <- seq_stats(Qy > high_thr)
      lp <- seq_stats(Qy < low_thr)
      
      # zero-flow percentile logic
      zp <- seq_stats(Qy <= zero_thr)
      
      # monthly means
      mtab <- tapply(Qy, mony, mean, na.rm = TRUE)
      mvec <- rep(NA_real_, 12)
      mvec[as.numeric(names(mtab))] <- mtab
      
      # rates
      dQy <- dQ[idx[-1]]  # diff eats first index
      
      res[[i]] <- c(
        amax1d   = max(Qy,   na.rm=TRUE),
        amax3d   = max(r3y,  na.rm=TRUE),
        amax7d   = max(r7y,  na.rm=TRUE),
        amax30d  = max(r30y, na.rm=TRUE),
        amax90d  = max(r90y, na.rm=TRUE),
        amin1d   = min(Qy,   na.rm=TRUE),
        amin3d   = min(r3y,  na.rm=TRUE),
        amin7d   = min(r7y,  na.rm=TRUE),
        amin30d  = min(r30y, na.rm=TRUE),
        amin90d  = min(r90y, na.rm=TRUE),
        
        juldate_amax1d = jdy[ which.max(Qy) ],
        juldate_amin1d = jdy[ which.min(Qy) ],
        
        num_high_pulses    = hp$num,
        avg_dur_high_pulses= hp$avg,
        num_low_pulses     = lp$num,
        avg_dur_low_pulses = lp$avg,
        
        mean_rate_rise = mean(dQy[dQy > 0],  na.rm=TRUE),
        mean_rate_fall = mean(abs(dQy[dQy < 0]), na.rm=TRUE),
        num_rises      = sum(dQy > 0),
        num_falls      = sum(dQy < 0),
        
        num_zero_seq    = zp$num,
        avg_dur_zero_seq= zp$avg,
        time_main_zero_seq =
          if (is.na(zp$start)) NA_real_ else jdy[zp$start],
        
        # monthly fields
        january_flow   = mvec[1],
        february_flow  = mvec[2],
        march_flow     = mvec[3],
        april_flow     = mvec[4],
        may_flow       = mvec[5],
        june_flow      = mvec[6],
        july_flow      = mvec[7],
        august_flow    = mvec[8],
        september_flow = mvec[9],
        october_flow   = mvec[10],
        november_flow  = mvec[11],
        december_flow  = mvec[12]
      )
    }
    
    # convert list → matrix
    mat <- do.call(rbind, res)
    
    # summary stats
    med <- apply(mat, 2, median, na.rm=TRUE)
    iqr <- apply(mat, 2, IQR,    na.rm=TRUE)
    
    out <- c(med, iqr)
    names(out) <- c(
      paste0(names(med), "_m"),
      paste0(names(iqr), "_v")
    )
    
    # metric metadata
    meta <- tibble::tribble(
      ~metric_id, ~group, ~metric,
      1,"HF","amax1d_m",2,"HF","amax1d_v",3,"HF","amax3d_m",4,"HF","amax3d_v",
      5,"HF","amax7d_m",6,"HF","amax7d_v",7,"HF","amax30d_m",8,"HF","amax30d_v",
      9,"HF","amax90d_m",10,"HF","amax90d_v",11,"HF","juldate_amax1d_m",12,"HF","juldate_amax1d_v",
      13,"HF","num_high_pulses_m",14,"HF","num_high_pulses_v",
      15,"HF","avg_dur_high_pulses_m",16,"HF","avg_dur_high_pulses_v",
      17,"MF","january_flow_m",18,"MF","january_flow_v",19,"MF","february_flow_m",
      20,"MF","february_flow_v",21,"MF","march_flow_m",22,"MF","march_flow_v",
      23,"MF","april_flow_m",24,"MF","april_flow_v",25,"MF","may_flow_m",
      26,"MF","may_flow_v",27,"MF","june_flow_m",28,"MF","june_flow_v",
      29,"MF","july_flow_m",30,"MF","july_flow_v",31,"MF","august_flow_m",
      32,"MF","august_flow_v",33,"MF","september_flow_m",
      34,"MF","september_flow_v",35,"MF","october_flow_m",36,"MF","october_flow_v",
      37,"MF","november_flow_m",38,"MF","november_flow_v",
      39,"MF","december_flow_m",40,"MF","december_flow_v",
      41,"LF","amin1d_m",42,"LF","amin1d_v",43,"LF","amin3d_m",44,"LF","amin3d_v",
      45,"LF","amin7d_m",46,"LF","amin7d_v",47,"LF","amin30d_m",48,"LF","amin30d_v",
      49,"LF","amin90d_m",50,"LF","amin90d_v",51,"LF","juldate_amin1d_m",
      52,"LF","juldate_amin1d_v",53,"LF","num_low_pulses_m",
      54,"LF","num_low_pulses_v",55,"LF","avg_dur_low_pulses_m",
      56,"LF","avg_dur_low_pulses_v",
      57,"RFC","mean_rate_rise_m",58,"RFC","mean_rate_rise_v",
      59,"RFC","mean_rate_fall_m",60,"RFC","mean_rate_fall_v",
      61,"RFC","num_rises_m",62,"RFC","num_rises_v",
      63,"RFC","num_falls_m",64,"RFC","num_falls_v",
      65,"IF","num_zero_seq_m",66,"IF","num_zero_seq_v",
      67,"IF","avg_dur_zero_seq_m",68,"IF","avg_dur_zero_seq_v",
      69,"IF","time_main_zero_seq_m",70,"IF","time_main_zero_seq_v"
    )
    
    tibble(
      metric_id = meta$metric_id,
      group     = factor(meta$group,
                         levels=c("HF","MF","LF","RFC","IF"), ordered=TRUE),
      metric    = meta$metric,
      value     = as.numeric(out[ meta$metric ])
    )
  }
  
  # --------------------------------------------------------------------------
  # MULTI-SITE HANDLING
  # --------------------------------------------------------------------------
  if (is.vector(Q)) {
    return(calc_site(Q, date))
  }
  
  if (!is.matrix(Q) && !is.data.frame(Q))
    stop("Q must be a vector, matrix, or data.frame")
  
  if (nrow(Q) != length(date))
    stop("nrow(Q) must equal length(date)")
  
  site_names <- colnames(Q)
  if (is.null(site_names))
    site_names <- paste0("site", seq_len(ncol(Q)))
  
  out_list <- vector("list", ncol(Q))
  
  for (i in seq_len(ncol(Q))) {
    site_out <- calc_site(Q[,i], date)
    site_out$site <- site_names[i]
    out_list[[i]] <- site_out
  }
  
  bind_rows(out_list) %>% relocate(site, .before=metric_id)
}
