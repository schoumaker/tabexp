#' Compute events and exposure  by years before survey (YBS)
#' from DHS-like data using person-period file logic.
#'
#' @import data.table
#' @importFrom haven read_dta zap_label zap_labels
#' @param data A data.frame or data.table with DHS variables.
#' @param ybs_max Maximum years before survey to include.
#' @param age_bins Age group breaks.
#' @param birth_vars Names of birth date variables (e.g., b3_01 to b3_20).
#' @param default_vars Names of essential variables (e.g., v011, v008, v005).
#'
#' @return A data.table with YBS, age group, events, exposure, and mid-period year.
#' @export



tabexp <- function(data = NULL,
                   ybs_max = 40,
                   age_bins = seq(10, 50, by = 5),
                   birth_vars = sprintf("b3_%02d", 1:20),
                   default_vars = c("caseid", "v011", "v008", "v005", "awfactt")) {

  library(data.table)
  library(haven)

  birth_vars <- sprintf("b3_%02d", 1:20)
  default_vars <- c("caseid", "v011", "v008", "v005", "awfactt")
  vars_to_load <- unique(c(default_vars, birth_vars))

  # --- Load data --
  raw <- as.data.table(data)

  raw[, weight := as.numeric(haven::zap_label(v005)) / 1e6]

  # Ensure awfactt exists
  if (!"awfactt" %in% names(raw)) {
    warning("`awfactt` not found in data. Defaulting to 100 (i.e., no adjustment for subsampling).")
    raw[, awfactt := 100]
  } else {
    message("Using `awfactt` to adjust for subsampling (weighted exposure = weight * awfactt / 100).")
  }


  raw[, awfactt := as.numeric(haven::zap_label(awfactt))]

  split_periods <- raw[, {
    age_at_interview <- v008 - v011
    max_months <- min(age_at_interview, ybs_max * 12)

    cutpoints <- sort(unique(c(
      v008 - seq(12, max_months, by = 12),
      v011 + seq(0, age_at_interview, by = 60)
    )))
    cutpoints <- cutpoints[cutpoints >= v011 & cutpoints < v008]
    start_cmc <- cutpoints
    end_cmc <- c(cutpoints[-1] - 1, v008 - 1)

    .(start_cmc = start_cmc,
      end_cmc = end_cmc,
      age = floor((start_cmc - v011) / 12),
      YBS = floor((v008 - start_cmc - 1) / 12),
      age_group = cut((start_cmc - v011) / 12, breaks = age_bins, right = FALSE),
      caseid = caseid,
      v005 = v005,
      weight = weight,
      v011 = v011,
      v008 = v008,
      awfactt = awfactt)
  }, by = caseid]

  births_long <- melt(raw,
                      id.vars = "caseid",
                      measure.vars = patterns("^b3_"),
                      value.name = "birth_cmc",
                      na.rm = TRUE)
  births_long[, `:=`(start = birth_cmc, end = birth_cmc)]

  intervals <- split_periods[, .(caseid, start = start_cmc, end = end_cmc, age_group, YBS)]
  setkey(intervals, caseid, start, end)

  events <- foverlaps(births_long, intervals, by.x = c("caseid", "start", "end"), nomatch = 0)


  birth_counts <- events[, .(births = .N), by = .(caseid, start, age_group, YBS)]
  setnames(birth_counts, "start", "start_cmc")

  split_periods <- split_periods[, .SD, .SDcols = !duplicated(names(split_periods))]
  split_periods <- merge(split_periods, birth_counts,
                         by = c("caseid", "start_cmc", "age_group", "YBS"),
                         all.x = TRUE)
  split_periods[is.na(births), births := 0]

  split_periods[, duration := (end_cmc - start_cmc + 1) / 12]
  split_periods[is.na(awfactt), awfactt := 100]
  split_periods[, py := duration * weight * awfactt / 100]
  split_periods[, events := births * weight]

  split_periods <- split_periods[!is.na(age_group)]

  agg_table_ybs_age <- split_periods[
    , .(
      births = sum(births, na.rm = TRUE),
      events = sum(events, na.rm = TRUE),
      exposure_py = sum(py, na.rm = TRUE)
    ),
    by = .(YBS, age_group)
  ][order(YBS, age_group)]

  # Convert YBS to mid-year using v008
  # Here we take the average interview date if multiple women
  mean_v008 <- mean(raw$v008, na.rm = TRUE)

  agg_table_ybs_age[, mid_cmc := mean_v008 - (YBS * 12 + 6)]
  agg_table_ybs_age[, mid_year := round(1900 + mid_cmc / 12, 1)]



  agg_table_ybs_age[, ASFR := events / exposure_py]

  result <- agg_table_ybs_age[, .(YBS, age_group, events, exposure_py, mid_year)]
  result <- haven::zap_labels(result)
  return(result)

}

