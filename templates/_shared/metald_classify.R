
count_cardiometabolic_risk <- function(bmi, waist, sex,
                                        fbg, hba1c, t2d_yes, t2d_med_yes,
                                        sbp, dbp, htn_med_yes,
                                        tg, lipid_med_yes,
                                        hdl) {
  c1 <- as.integer(
    (!is.na(bmi) & bmi >= 25) |
    (!is.na(waist) & ((sex == "Male" & waist > 94) | (sex == "Female" & waist > 80)))
  )
  c2 <- as.integer(
    (!is.na(fbg) & fbg >= 100) |
    (!is.na(hba1c) & hba1c >= 5.7) |
    (!is.na(t2d_yes) & t2d_yes == 1) |
    (!is.na(t2d_med_yes) & t2d_med_yes == 1)
  )
  c3 <- as.integer(
    (!is.na(sbp) & sbp >= 130) |
    (!is.na(dbp) & dbp >= 85) |
    (!is.na(htn_med_yes) & htn_med_yes == 1)
  )
  c4 <- as.integer(
    (!is.na(tg) & tg >= 150) |
    (!is.na(lipid_med_yes) & lipid_med_yes == 1)
  )
  c5 <- as.integer(
    (!is.na(hdl) & ((sex == "Male" & hdl < 40) | (sex == "Female" & hdl < 50))) |
    (!is.na(lipid_med_yes) & lipid_med_yes == 1)
  )
  rowSums(cbind(c1, c2, c3, c4, c5), na.rm = TRUE)
}

alcohol_gwk_from_alq130 <- function(alq130) {
  ifelse(is.na(alq130) | alq130 %in% c(77, 99, 777, 999), NA, alq130 * 14 * 7)
}

fli_threelevel <- function(fli) {
  factor(dplyr::case_when(
    is.na(fli) ~ NA_character_,
    fli < 30   ~ "No steatosis (<30)",
    fli < 60   ~ "Indeterminate (30-60)",
    fli >= 60  ~ "Steatosis (>=60)"
  ), levels = c("No steatosis (<30)", "Indeterminate (30-60)", "Steatosis (>=60)"))
}

metald_classify <- function(steatosis_yes, alcohol_gwk, sex, cm_risk_count) {

  thresh_low  <- ifelse(sex == "Female", 140, 210)
  thresh_high <- ifelse(sex == "Female", 350, 420)

  factor(
    dplyr::case_when(
      is.na(steatosis_yes)                  ~ NA_character_,
      steatosis_yes == 0                    ~ "No steatosis",
      steatosis_yes == 1 & cm_risk_count == 0 & (!is.na(alcohol_gwk) & alcohol_gwk > thresh_high)
                                            ~ "ALD",
      steatosis_yes == 1 & cm_risk_count >= 1 & (is.na(alcohol_gwk) | alcohol_gwk < thresh_low)
                                            ~ "MASLD",
      steatosis_yes == 1 & cm_risk_count >= 1 & alcohol_gwk >= thresh_low & alcohol_gwk <= thresh_high
                                            ~ "MetALD",
      steatosis_yes == 1 & cm_risk_count >= 1 & alcohol_gwk > thresh_high
                                            ~ "ALD",
      steatosis_yes == 1 & cm_risk_count == 0
                                            ~ "Cryptogenic steatosis",
      TRUE                                  ~ NA_character_
    ),
    levels = c("No steatosis", "MASLD", "MetALD", "ALD", "Cryptogenic steatosis")
  )
}
