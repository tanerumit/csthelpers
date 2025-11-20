#' Calculate Environmental Flow (E-flow) Hydrologic Metrics
#'
#' @description
#' Computes 70 eco-hydrological indicators describing magnitude, timing,
#' frequency, rate of change, and intermittency/zero-flow characteristics of
#' discharge time-series. Supports both single-site and multi-site input.
#'
#' @param Q Numeric vector or matrix/data.frame of daily discharge values.
#'   If a vector, a single site is processed. If a matrix/data.frame, each
#'   column is treated as a separate site.
#' @param date Vector of corresponding dates (`Date` or coercible). Length
#'   must match the number of rows in `Q`.
#' @param high_thresh_quantile Quantile used to define high-flow pulse
#'   threshold. Default: `0.9`.
#' @param low_thresh_quantile Quantile used to define low-flow pulse
#'   threshold. Default: `0.1`.
#' @param zero_thresh_quantile Quantile defining the "near-zero" threshold
#'   used for IF (intermittency) metrics. Default: `0.02`.
#' 
#' @details
#' For each site, the function:
#' \itemize{
#'   \item Derives 35 yearly metrics (magnitude, timing, pulses,
#'         rate of change, zero-flow).
#'   \item Aggregates them to median (`*_m`) and IQR (`*_v`) across years.
#'   \item Returns a tidy table with 70 indicators grouped into HF, MF, LF,
#'         RFC, IF.
#' }
#'
#' IF-group metrics (zero-flow count, duration, timing) use flows
#' \eqn{Q \le \text{quantile}(Q, zero\_thresh\_quantile)} rather than Q == 0.
#'
#' @return A tibble containing:
#'   \itemize{
#'     \item site (if multi-site)
#'     \item metric_id
#'     \item group (HF, MF, LF, RFC, IF)
#'     \item metric (short indicator name)
#'     \item value (numeric)
#'   }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' d <- seq.Date(as.Date("1990-01-01"), as.Date("2010-12-31"), by="day")
#' Q <- pmax(0, 30 + 10*sin(2*pi*(1:length(d))/365) + rnorm(length(d), 0, 3))
#'
#' res <- calculate_eflow_metrics(Q, d, zero_thresh_quantile = 0.02)
#' head(res)
#' }
#'
#' @import dplyr lubridate zoo tibble
#' @export
calculate_eflow_metrics_old <- function(
    Q, date,
    high_thresh_quantile = 0.9,
    low_thresh_quantile  = 0.1,
    zero_thresh_quantile = 0.02
) {
  
  # -------------------------
  # Helper for single site
  # -------------------------
  calc_site <- function(Qv, datev) {
    
    stopifnot(length(Qv) == length(datev))
    if (!inherits(datev, "Date")) datev <- as.Date(datev)
    
    df <- data.frame(date = datev, Q = as.numeric(Qv))
    
    # Internal helper functions
    rollmean_n <- function(x, n) zoo::rollmean(x, k = n, align = "right", fill = NA)
    
    seq_stats <- function(cond) {
      r <- rle(cond)
      n <- sum(r$values)
      dur <- if (n > 0) mean(r$lengths[r$values]) else 0
      list(num = n, avg_dur = dur)
    }
    
    # Preprocessing
    df <- df %>%
      mutate(
        year = year(date),
        month = month(date),
        jday = yday(date)
      ) %>%
      arrange(date)
    
    high_thr <- quantile(df$Q, probs = high_thresh_quantile, na.rm = TRUE)
    low_thr  <- quantile(df$Q, probs = low_thresh_quantile,  na.rm = TRUE)
    zero_thr <- quantile(df$Q, probs = zero_thresh_quantile, na.rm = TRUE)
    
    yrs <- unique(df$year)
    res <- data.frame(year = yrs)
    
    # Loop through years
    for (y in yrs) {
      sub <- df[df$year == y, ]
      if (nrow(sub) < 30) next
      
      sub$roll3  <- rollmean_n(sub$Q, 3)
      sub$roll7  <- rollmean_n(sub$Q, 7)
      sub$roll30 <- rollmean_n(sub$Q, 30)
      sub$roll90 <- rollmean_n(sub$Q, 90)
      
      # Magnitude
      res[res$year == y, "amax1d"]  <- max(sub$Q, na.rm = TRUE)
      res[res$year == y, "amax3d"]  <- max(sub$roll3, na.rm = TRUE)
      res[res$year == y, "amax7d"]  <- max(sub$roll7, na.rm = TRUE)
      res[res$year == y, "amax30d"] <- max(sub$roll30, na.rm = TRUE)
      res[res$year == y, "amax90d"] <- max(sub$roll90, na.rm = TRUE)
      
      res[res$year == y, "amin1d"]  <- min(sub$Q, na.rm = TRUE)
      res[res$year == y, "amin3d"]  <- min(sub$roll3, na.rm = TRUE)
      res[res$year == y, "amin7d"]  <- min(sub$roll7, na.rm = TRUE)
      res[res$year == y, "amin30d"] <- min(sub$roll30, na.rm = TRUE)
      res[res$year == y, "amin90d"] <- min(sub$roll90, na.rm = TRUE)
      
      # Timing
      res[res$year == y, "juldate_amax1d"] <- sub$jday[which.max(sub$Q)]
      res[res$year == y, "juldate_amin1d"] <- sub$jday[which.min(sub$Q)]
      
      # High/Low pulses
      h <- seq_stats(sub$Q > high_thr)
      l <- seq_stats(sub$Q < low_thr)
      
      res[res$year == y, c("num_high_pulses","avg_dur_high_pulses")] <- c(h$num, h$avg_dur)
      res[res$year == y, c("num_low_pulses","avg_dur_low_pulses")]   <- c(l$num, l$avg_dur)
      
      # Rate of change
      dQ <- diff(sub$Q)
      res[res$year == y, "mean_rate_rise"] <- mean(dQ[dQ > 0], na.rm = TRUE)
      res[res$year == y, "mean_rate_fall"] <- mean(abs(dQ[dQ < 0]), na.rm = TRUE)
      res[res$year == y, "num_rises"] <- sum(dQ > 0, na.rm = TRUE)
      res[res$year == y, "num_falls"] <- sum(dQ < 0, na.rm = TRUE)
      
      # -------------------------------------------------
      # IF metrics using threshold: Q <= zero_thr
      # -------------------------------------------------
      is_zero <- sub$Q <= zero_thr
      
      z <- seq_stats(is_zero)
      res[res$year == y, c("num_zero_seq","avg_dur_zero_seq")] <- c(z$num, z$avg_dur)
      
      if (z$num > 0) {
        r <- rle(is_zero)
        zi <- which(r$values)
        main_i <- zi[which.max(r$lengths[zi])]
        start <- sum(r$lengths[1:(main_i - 1)]) + 1
        res[res$year == y, "time_main_zero_seq"] <- sub$jday[start]
      } else {
        res[res$year == y, "time_main_zero_seq"] <- NA
      }
      
      # Monthly means
      mtab <- sub %>%
        group_by(month) %>%
        summarise(flow = mean(Q, na.rm = TRUE), .groups = "drop")
      
      for (m in 1:12) {
        nm <- paste0(tolower(month.name[m]), "_flow")
        res[res$year == y, nm] <- if (m %in% mtab$month) mtab$flow[mtab$month == m] else NA
      }
    }
    
    # Summaries
    stats <- res %>% select(-year)
    med <- apply(stats, 2, median, na.rm = TRUE)
    iqr <- apply(stats, 2, IQR,    na.rm = TRUE)
    
    out <- as.data.frame(as.list(c(
      setNames(med, paste0(names(med), "_m")),
      setNames(iqr, paste0(names(iqr), "_v"))
    )))
    
    # Classification / metadata table
    group_tbl <- tibble::tribble(
      ~metric_id, ~group, ~metric,
      1,"HF","amax1d_m",2,"HF","amax1d_v",3,"HF","amax3d_m",4,"HF","amax3d_v",
      5,"HF","amax7d_m",6,"HF","amax7d_v",7,"HF","amax30d_m",8,"HF","amax30d_v",
      9,"HF","amax90d_m",10,"HF","amax90d_v",11,"HF","juldate_amax1d_m",12,"HF","juldate_amax1d_v",
      13,"HF","num_high_pulses_m",14,"HF","num_high_pulses_v",15,"HF","avg_dur_high_pulses_m",
      16,"HF","avg_dur_high_pulses_v",
      17,"MF","january_flow_m",18,"MF","january_flow_v",19,"MF","february_flow_m",
      20,"MF","february_flow_v",21,"MF","march_flow_m",22,"MF","march_flow_v",
      23,"MF","april_flow_m",24,"MF","april_flow_v",25,"MF","may_flow_m",26,"MF","may_flow_v",
      27,"MF","june_flow_m",28,"MF","june_flow_v",29,"MF","july_flow_m",30,"MF","july_flow_v",
      31,"MF","august_flow_m",32,"MF","august_flow_v",33,"MF","september_flow_m",34,"MF","september_flow_v",
      35,"MF","october_flow_m",36,"MF","october_flow_v",37,"MF","november_flow_m",38,"MF","november_flow_v",
      39,"MF","december_flow_m",40,"MF","december_flow_v",
      41,"LF","amin1d_m",42,"LF","amin1d_v",43,"LF","amin3d_m",44,"LF","amin3d_v",
      45,"LF","amin7d_m",46,"LF","amin7d_v",47,"LF","amin30d_m",48,"LF","amin30d_v",
      49,"LF","amin90d_m",50,"LF","amin90d_v",51,"LF","juldate_amin1d_m",52,"LF","juldate_amin1d_v",
      53,"LF","num_low_pulses_m",54,"LF","num_low_pulses_v",55,"LF","avg_dur_low_pulses_m",56,"LF","avg_dur_low_pulses_v",
      57,"RFC","mean_rate_rise_m",58,"RFC","mean_rate_rise_v",59,"RFC","mean_rate_fall_m",60,"RFC","mean_rate_fall_v",
      61,"RFC","num_rises_m",62,"RFC","num_rises_v",63,"RFC","num_falls_m",64,"RFC","num_falls_v",
      65,"IF","num_zero_seq_m",66,"IF","num_zero_seq_v",67,"IF","avg_dur_zero_seq_m",68,"IF","avg_dur_zero_seq_v",
      69,"IF","time_main_zero_seq_m",70,"IF","time_main_zero_seq_v"
    )
    
    tidy_out <- group_tbl %>%
      mutate(value = as.numeric(out[1, metric])) %>%
      mutate(group = factor(group,
                            levels = c("HF","MF","LF","RFC","IF"),
                            ordered = TRUE))
    
    return(tidy_out)
  }
  
  # -------------------------
  # Multi-site handling
  # -------------------------
  if (is.vector(Q)) {
    
    res <- calc_site(Q, date)
    
  } else if (is.matrix(Q) || is.data.frame(Q)) {
    
    if (nrow(Q) != length(date))
      stop("If Q is a matrix/data.frame, number of rows must equal length(date).")
    
    site_names <- if (!is.null(colnames(Q))) colnames(Q) else paste0("site", seq_len(ncol(Q)))
    
    site_results <- lapply(seq_len(ncol(Q)), function(i) {
      tbl <- calc_site(Q[, i], date)
      tbl$site <- site_names[i]
      tbl
    })
    
    res <- bind_rows(site_results) %>%
      relocate(site, .before = metric_id)
    
  } else {
    stop("Q must be a numeric vector, matrix, or data frame.")
  }
  
  return(res)
}
