
suppressPackageStartupMessages({
  library(dplyr)
})


nhanes_cycle_table <- function() {
  data.frame(
    cycle_tag = c(
      "NHANES_1999_2000", "NHANES_2001_2002", "NHANES_2003_2004",
      "NHANES_2005_2006", "NHANES_2007_2008", "NHANES_2009_2010",
      "NHANES_2011_2012", "NHANES_2013_2014", "NHANES_2015_2016",
      "NHANES_2017_2018", "PrePandemic_2017_March2020"
    ),
    suffix = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "P"),
    pathYear = c(1999L, 2001L, 2003L, 2005L, 2007L, 2009L,
                 2011L, 2013L, 2015L, 2017L, 2017L),
    years_label = c(
      "1999-2000", "2001-2002", "2003-2004", "2005-2006",
      "2007-2008", "2009-2010", "2011-2012", "2013-2014",
      "2015-2016", "2017-2018", "2017-March 2020"
    ),
    years = c(1999L, 2001L, 2003L, 2005L, 2007L, 2009L,
              2011L, 2013L, 2015L, 2017L, 2017L),
    span_yr = c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3.2),
    weight_mec   = c(rep("WTMEC2YR", 10), "WTMECPRP"),
    weight_diet1 = c(rep("WTDRD1", 10), "WTDRD1PP"),
    weight_diet2 = c(rep("WTDR2D", 10), "WTDR2DPP"),
    weight_saf   = c(rep("WTSAF2YR", 10), "WTSAFPRP"),
    weight_int   = c(rep("WTINT2YR", 10), "WTINTPRP"),
    is_prepandemic = c(rep(FALSE, 10), TRUE),
    stringsAsFactors = FALSE
  )
}


nhanes_modules_standard <- function() {
  list(
    core = c("DEMO", "BMX", "BPX", "ALQ", "DBQ", "SMQ",
             "PAQ", "RHQ", "DIQ", "BPQ", "MCQ"),
    metals = c("PBCD", "IHGEM", "UTAS", "UAS", "UHM"),
    fibroscan = c("LUX"),
    diet = c("DR1TOT", "DR2TOT"),
    biochem = c("BIOPRO", "GHB", "GLU", "TCHOL", "HDL",
                "TRIGLY", "INS", "CBC", "ALB_CR",
                "HEPB_S", "HEPC"),
    mortality = c("MORT"),
    rx = c("RXQ_RX")
  )
}


make_cycle_spec <- function(year_range,
                            include_prepandemic = FALSE,
                            exclude_cycles = character(0)) {
  ct <- nhanes_cycle_table()

  ct_std <- subset(ct, !is_prepandemic & years %in% year_range)

  ct_pp <- if (include_prepandemic) subset(ct, is_prepandemic) else ct[FALSE, ]

  ct_keep <- rbind(ct_std, ct_pp)
  ct_keep <- subset(ct_keep, !(cycle_tag %in% exclude_cycles))

  if (nrow(ct_keep) == 0) {
    stop("make_cycle_spec: no cycles match year_range=", paste(range(year_range), collapse = ":"),
         ", include_prepandemic=", include_prepandemic)
  }

  out <- lapply(seq_len(nrow(ct_keep)), function(i) {
    list(
      suffix         = ct_keep$suffix[i],
      pathYear       = ct_keep$pathYear[i],
      years_label    = ct_keep$years_label[i],
      span_yr        = ct_keep$span_yr[i],
      weight_mec     = ct_keep$weight_mec[i],
      weight_diet1   = ct_keep$weight_diet1[i],
      weight_diet2   = ct_keep$weight_diet2[i],
      weight_saf     = ct_keep$weight_saf[i],
      weight_int     = ct_keep$weight_int[i],
      is_prepandemic = ct_keep$is_prepandemic[i]
    )
  })
  names(out) <- ct_keep$cycle_tag
  out
}


cycle_years_vector <- function(cycles_list) {
  yrs <- vapply(cycles_list, function(x) x$span_yr, numeric(1))
  names(yrs) <- names(cycles_list)
  yrs
}


validate_no_J_in_PP <- function(cycles_list) {
  has_J  <- "NHANES_2017_2018" %in% names(cycles_list)
  has_PP <- "PrePandemic_2017_March2020" %in% names(cycles_list)
  if (has_J && has_PP) {
    warning("[cycle_config] Both NHANES_2017_2018 (J) and PrePandemic_2017_March2020 (P_) ",
            "present in cycle spec. J ⊂ P_ — you must dedupe (drop J rows where SEQN is ",
            "also in P_, OR use P_ only). See NCHS Series 2 No. 190 §4. Examples: ",
            "main = P_ only; main = G+H+I+J no P_.")
    return(invisible(FALSE))
  }
  invisible(TRUE)
}