
library(dplyr)


COV_BASIC <- c("age", "sex", "race", "education", "marital", "pir")

COV_M1 <- c(COV_BASIC, "smoke", "drink")

COV_M2_COX <- c(COV_M1, "bmi", "fbg", "hba1c", "ldl",
                "sbp", "dbp", "diabetes", "hypertension",
                "family_hx_heart")

COV_M2_LOGIT <- c(COV_M2_COX, "alt", "alp", "albumin")

COV_PREGNANCY <- c("age", "race", "education", "marital", "pir",
                   "pre_pregnancy_bmi", "smoke", "drink",
                   "physical_activity", "vitamin_d", "parity")

SUBGROUP_VARS <- c("sex", "age_group", "race", "education",
                   "diabetes", "hypertension", "pir_group",
                   "smoke", "drink")


education_recode <- function(x) {
  case_when(
    x %in% c(1, 2) ~ "Less than high school",
    x == 3 ~ "High school",
    x %in% c(4, 5) ~ "College or above",
    TRUE ~ NA_character_
  )
}

marital_recode <- function(x) {
  case_when(
    x %in% c(1, 6) ~ "Married/Living with partner",
    x %in% c(2, 3, 4, 5) ~ "Not married",
    TRUE ~ NA_character_
  )
}

pir_recode <- function(x) {
  case_when(
    x <= 1.3 ~ "Low",
    x <= 3.5 ~ "Middle",
    x > 3.5 ~ "High",
    TRUE ~ NA_character_
  )
}

smoke_recode <- function(smq020, smq040) {
  case_when(
    smq020 == 2 ~ "Never",
    smq020 == 1 & smq040 == 3 ~ "Former",
    smq020 == 1 & smq040 %in% c(1, 2) ~ "Current",
    TRUE ~ NA_character_
  )
}

diabetes_define <- function(diq, diq050, diq070, lbxglu, lbxgh) {
  case_when(
    diq== 1 | diq050 == 1 | diq070 == 1 ~ "Yes",
    !is.na(lbxglu) & lbxglu >= 7 ~ "Yes",
    !is.na(lbxgh) & lbxgh >= 6.5 ~ "Yes",
    TRUE ~ "No"
  )
}

hypertension_define <- function(bpq020, sbp_mean, dbp_mean) {
  case_when(
    bpq020 == 1 | sbp_mean >= 140 | dbp_mean >= 90 ~ "Yes",
    TRUE ~ "No"
  )
}

coalesce_cols <- function(df, cols) {
  cols <- intersect(cols, names(df))
  if (length(cols) == 0) return(rep(NA_real_, nrow(df)))
  if (length(cols) == 1) return(df[[cols[1]]])
  Reduce(function(a, b) ifelse(is.na(a), b, a),
         lapply(cols, function(c) df[[c]]))
}

NHANES_VAR_FALLBACK <- list(
  hdl_mgdl     = c("LBDHDL", "LBXHDD", "LBDHDD"),
  sleep_hours  = c("SLDH", "SLD"),
  vit_d_serum  = c("LBXVIDMS", "LBXVDMS", "LBXVD2MS", "LBXVD3MS"),
  triglyceride = c("LBXTR"),
  alt          = c("LBXSATSI", "LBXSALT"),
  ast          = c("LBXSASSI", "LBXSAST"),
  hba1c        = c("LBXGH"),
  insulin      = c("LBXIN")
)
