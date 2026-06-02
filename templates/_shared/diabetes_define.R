
suppressPackageStartupMessages({
  library(dplyr)
})


ANTIDIABETIC_DRUG_PATTERNS <- c(
  "metformin", "glucophage", "fortamet", "riomet", "glumetza",
  "glipizide", "glyburide", "glibenclamide", "glimepiride",
  "tolazamide", "tolbutamide", "chlorpropamide", "amaryl", "glucotrol", "diabeta", "micronase",
  "repaglinide", "nateglinide", "prandin", "starlix",
  "pioglitazone", "rosiglitazone", "actos", "avandia",
  "sitagliptin", "saxagliptin", "linagliptin", "alogliptin",
  "januvia", "onglyza", "tradjenta", "nesina",
  "empagliflozin", "dapagliflozin", "canagliflozin", "ertugliflozin",
  "jardiance", "farxiga", "invokana", "steglatro",
  "liraglutide", "semaglutide", "dulaglutide", "exenatide", "lixisenatide", "tirzepatide",
  "victoza", "ozempic", "rybelsus", "trulicity", "byetta", "bydureon", "adlyxin", "mounjaro",
  "wegovy", "zepbound",
  "acarbose", "miglitol", "precose", "glyset",
  "insulin",
  "glargine", "detemir", "degludec", "lispro", "aspart", "glulisine",
  "lantus", "basaglar", "toujeo", "levemir", "tresiba",
  "humalog", "novolog", "apidra", "humulin", "novolin",
  "pramlintide", "symlin",
  "janumet", "kombiglyze", "jentadueto", "kazano", "synjardy",
  "xigduo", "invokamet", "segluromet", "soliqua", "xultophy"
)


detect_antidiabetic_rx <- function(rxq_rx_df, drug_patterns = ANTIDIABETIC_DRUG_PATTERNS) {
  if (!all(c("SEQN", "RXDDRUG") %in% names(rxq_rx_df))) {
    warning("[detect_antidiabetic_rx] rxq_rx_df missing SEQN or RXDDRUG; returning empty")
    return(data.frame(SEQN = integer(0), rx_diabetes_yes = logical(0)))
  }
  drug_lc <- tolower(as.character(rxq_rx_df$RXDDRUG))
  pat <- paste(drug_patterns, collapse = "|")
  hit <- grepl(pat, drug_lc, perl = TRUE)
  agg <- aggregate(hit, by = list(SEQN = rxq_rx_df$SEQN), FUN = function(x) any(x, na.rm = TRUE))
  names(agg) <- c("SEQN", "rx_diabetes_yes")
  agg
}


diabetes_comprehensive <- function(df) {

  get_col <- function(canonical, raw, df) {
    if (canonical %in% names(df)) return(df[[canonical]])
    if (raw       %in% names(df)) return(df[[raw]])
    return(rep(NA_real_, nrow(df)))
  }

  diq   <- get_col("diabetes_self",    "DIQ", df)
  diq050    <- get_col("diabetes_insulin", "DIQ050", df)
  diq070    <- get_col("diabetes_oralmed", "DIQ070", df)
  hba1c     <- get_col("hba1c_pct",        "LBXGH",  df)
  fpg       <- get_col("fpg_mgdl",         "LBXGLU", df)
  ogtt_2h   <- get_col("ogtt_2hpg_mgdl",   "LBXGLT", df)
  rx_hit <- if ("rx_diabetes_yes" %in% names(df)) {
    as.logical(df$rx_diabetes_yes)
  } else {
    rep(NA, nrow(df))
  }

  ev_self    <- !is.na(diq) & diq== 1
  ev_insulin <- !is.na(diq050) & diq050 == 1
  ev_oral    <- !is.na(diq070) & diq070 == 1
  ev_a1c     <- !is.na(hba1c)  & hba1c  >= 6.5
  ev_fpg     <- !is.na(fpg)    & fpg    >= 126
  ev_ogtt    <- !is.na(ogtt_2h) & ogtt_2h >= 200
  ev_rx      <- !is.na(rx_hit) & rx_hit

  is_dm <- ev_self | ev_insulin | ev_oral | ev_a1c | ev_fpg | ev_ogtt | ev_rx

  has_any_signal <- !is.na(diq) | !is.na(diq050) | !is.na(diq070) |
                    !is.na(hba1c)  | !is.na(fpg)    | !is.na(ogtt_2h) |
                    !is.na(rx_hit)

  out <- rep(NA_integer_, nrow(df))
  out[is_dm] <- 1L
  out[!is_dm & has_any_signal] <- 0L

  out
}


diabetes_criteria_breakdown <- function(df) {
  get_col <- function(canonical, raw, df) {
    if (canonical %in% names(df)) return(df[[canonical]])
    if (raw       %in% names(df)) return(df[[raw]])
    return(rep(NA_real_, nrow(df)))
  }

  diq  <- get_col("diabetes_self",    "DIQ", df)
  diq050   <- get_col("diabetes_insulin", "DIQ050", df)
  diq070   <- get_col("diabetes_oralmed", "DIQ070", df)
  hba1c    <- get_col("hba1c_pct",        "LBXGH",  df)
  fpg      <- get_col("fpg_mgdl",         "LBXGLU", df)
  ogtt_2h  <- get_col("ogtt_2hpg_mgdl",   "LBXGLT", df)
  rx_hit <- if ("rx_diabetes_yes" %in% names(df)) df$rx_diabetes_yes else rep(NA, nrow(df))

  n_total <- nrow(df)
  rows <- list(
    list("DIQ== 1 (self-report Yes)",       sum(!is.na(diq) & diq== 1)),
    list("DIQ050 == 1 (insulin)",                sum(!is.na(diq050) & diq050 == 1)),
    list("DIQ070 == 1 (oral hypoglycemic)",      sum(!is.na(diq070) & diq070 == 1)),
    list("LBXGH >= 6.5 (HbA1c)",                 sum(!is.na(hba1c)  & hba1c  >= 6.5)),
    list("LBXGLU >= 126 (fasting glucose)",      sum(!is.na(fpg)    & fpg    >= 126)),
    list("LBXGLT >= 200 (2h-OGTT)",              sum(!is.na(ogtt_2h) & ogtt_2h >= 200)),
    list("RXQ_RX antidiabetic drug",             sum(!is.na(rx_hit) & rx_hit))
  )
  do.call(rbind, lapply(rows, function(r) {
    data.frame(
      criterion  = r[[1]],
      n_positive = r[[2]],
      pct_of_total = round(100 * r[[2]] / n_total, 2),
      stringsAsFactors = FALSE
    )
  }))
}