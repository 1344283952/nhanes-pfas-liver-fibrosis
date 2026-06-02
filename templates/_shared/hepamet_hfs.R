
calc_hepamet_hfs <- function(sex_male, age, ast, albumin, platelet,
                              homa_ir, diabetes_yes) {

  z <- -5.390 +
       0.986 * sex_male +
       0.875 * as.integer(age >= 45) +
       0.012 * ast +
      -0.808 * albumin +                     # FIXED: continuous, not categorical
      -0.011 * platelet +
       0.072 * homa_ir +
       1.302 * diabetes_yes

  hfs <- 1 / (1 + exp(-z))
  hfs
}

hepamet_category <- function(hfs) {
  factor(
    cut(hfs, breaks = c(-Inf, 0.12, 0.47, Inf),
        labels = c("Low risk (<0.12)", "Indeterminate (0.12-0.47)", "High risk (>=0.47)"),
        right = FALSE, include.lowest = TRUE),
    levels = c("Low risk (<0.12)", "Indeterminate (0.12-0.47)", "High risk (>=0.47)")
  )
}

