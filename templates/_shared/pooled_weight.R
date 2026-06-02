
pooled_mec_weight <- function(df, cycle_years) {
  total_years <- sum(cycle_years)
  if (!"cycle_tag" %in% names(df)) {
    stop("df must have cycle_tag column for pooled_mec_weight")
  }
  yrs_per_row <- cycle_years[as.character(df$cycle_tag)]
  if (any(is.na(yrs_per_row))) {
    warning("Unknown cycle_tag values in df; weights set to NA")
  }
  is_prepandemic <- grepl("PrePandemic|^P_", as.character(df$cycle_tag))
  base_weight <- ifelse(is_prepandemic, df$WTMECPRP, df$WTMEC2YR)
  base_weight * (yrs_per_row / total_years)
}

pooled_saf_weight <- function(df, cycle_years) {
  total_years <- sum(cycle_years)
  yrs_per_row <- cycle_years[as.character(df$cycle_tag)]
  is_prepandemic <- grepl("PrePandemic|^P_", as.character(df$cycle_tag))
  base_weight <- ifelse(is_prepandemic,
                        df$WTSAFPRP %||% df$WTMECPRP,
                        df$WTSAF2YR %||% df$WTMEC2YR)
  base_weight * (yrs_per_row / total_years)
}

pooled_diet_weight <- function(df, cycle_years) {
  total_years <- sum(cycle_years)
  yrs_per_row <- cycle_years[as.character(df$cycle_tag)]
  is_prepandemic <- grepl("PrePandemic|^P_", as.character(df$cycle_tag))
  base_weight <- ifelse(is_prepandemic,
                        df$WTDRD1PP %||% df$WTMECPRP,
                        df$WTDRD1 %||% df$WTMEC2YR)
  base_weight * (yrs_per_row / total_years)
}

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
