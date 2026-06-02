
calc_fib4 <- function(age, ast, alt, plt) {
  fib4 <- (age * ast) / (plt * sqrt(alt))
  fib4
}

fib4_category <- function(fib4) {
  factor(
    cut(fib4, breaks = c(-Inf, 1.30, 2.67, Inf),
        labels = c("Low risk (<1.30)", "Indeterminate (1.30-2.67)", "High risk (>=2.67)"),
        right = FALSE, include.lowest = TRUE),
    levels = c("Low risk (<1.30)", "Indeterminate (1.30-2.67)", "High risk (>=2.67)")
  )
}

fib4_binary_advanced <- function(fib4) as.integer(!is.na(fib4) & fib4 >= 1.30)

fib4_advanced_age_adj <- function(fib4, age) {
  cut <- ifelse(age >= 65, 2.0, 1.30)
  as.integer(!is.na(fib4) & fib4 >= cut)
}

fib4_high_risk <- function(fib4) as.integer(!is.na(fib4) & fib4 >= 2.67)

