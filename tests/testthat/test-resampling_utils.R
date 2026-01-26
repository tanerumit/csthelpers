# Functions: find_nn1 (R/resampling_utils.R), select_maximin (R/resampling_utils.R), slice_event_window (R/resampling_utils.R), find_sequences_below (R/resampling_utils.R)

testthat::test_that("find_nn1 returns nearest indices and distances", {
  data <- data.frame(x = c(0, 1, 2), y = c(0, 1, 2))
  query <- data.frame(x = c(0.1, 1.8), y = c(0.1, 1.9))

  idx <- find_nn1(data, query)
  testthat::expect_equal(idx, c(1L, 3L))

  res <- find_nn1(data, query, return_dist = TRUE)
  testthat::expect_equal(res$idx, c(1L, 3L))
  testthat::expect_true(all(res$dist >= 0))
})

testthat::test_that("select_maximin returns k indices within range", {
  x <- c(0, 0, 1, 1)
  y <- c(0, 1, 0, 1)
  idx <- select_maximin(x, y, k = 2)
  testthat::expect_equal(length(idx), 2)
  testthat::expect_true(all(idx %in% seq_along(x)))
})

testthat::test_that("slice_event_window returns fixed-length windows", {
  events <- as.Date(c("2005-06-10", "2005-06-12"))
  timeline <- seq.Date(as.Date("2000-01-01"), as.Date("2010-12-31"), by = "day")

  idx <- slice_event_window(
    events = events,
    timeline = timeline,
    days.before = 365,
    total.days = 365,
    month.start = 1
  )

  testthat::expect_equal(length(idx), 365)
  testthat::expect_true(all(timeline[idx] >= as.Date("2000-01-01")))

  idx_list <- slice_event_window(
    events = list(events, as.Date("2008-02-01")),
    timeline = timeline,
    days.before = 200,
    total.days = 200,
    month.start = 6
  )
  testthat::expect_true(is.list(idx_list))
  testthat::expect_equal(length(idx_list[[1]]), length(idx_list[[2]]))
})

testthat::test_that("find_sequences_below finds runs and respects min_start_index", {
  x <- c(5, 2, 1, 1, 5, 2, 2, 2, 2)
  res <- find_sequences_below(x, threshold = 3, min_length = 3, min_start_index = 1)
  testthat::expect_equal(res, list(2:4, 6:9))

  res2 <- find_sequences_below(x, threshold = 3, min_length = 3, min_start_index = 7)
  testthat::expect_equal(res2, list())
})
