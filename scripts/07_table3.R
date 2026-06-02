# ============================================
# scripts/07_table3.R → qgcomp + WQS PFAS mixture triangulation
# Stack 1 主 (FIB-4 outcome) + Stack 7 (mortality)
# 加项 #7+8: Mixture exposure 多方法 triangulation
# ============================================

suppressPackageStartupMessages({
  library(qgcomp); library(gWQS); library(dplyr); library(survey)
})

cat("========================================\n")
cat("qgcomp + WQS PFAS mixture (#7+8)\n")
cat("========================================\n\n")

load("data/processed/nhanes_final.RData")
source("templates/_shared/covariates_config.R")
set.seed(20260519)

# 6 PFAS imputed
pfas_cols <- intersect(c("LBXPFOA_imp","LBXPFOS_imp","LBXPFNA_imp","LBXPFHS_imp",
                          "LBXPFDE_imp","LBXMPAH_imp"),
                       names(nhanes_final))
cat(sprintf("PFAS mixture: %d compounds (%s)\n",
            length(pfas_cols), paste(pfas_cols, collapse=", ")))

# Stack 1: outcome = fib4_log
# 用 minimal numeric/binary cov set (避免 factor cov subset 后 1-level + fasting-only NA)
cov_minimal <- intersect(c("age","RIAGENDR","race","education","pir","smoke","bmi","diabetes","hypertension"),
                         names(nhanes_final))
df1 <- nhanes_final %>%
  filter(if_all(all_of(pfas_cols), ~ !is.na(.)),
         !is.na(fib4_log)) %>%
  filter(if_all(all_of(cov_minimal), ~ !is.na(.))) %>%
  mutate(across(where(is.factor), droplevels))   # drop empty factor levels (fix WQS contrasts)
cat(sprintf("Stack 1 qgcomp+WQS cohort n=%d (minimal numeric cov set)\n", nrow(df1)))

cov_use <- cov_minimal
f_lhs <- "fib4_log"
f_rhs <- paste(c(pfas_cols, cov_use), collapse = " + ")

# qgcomp Stack 1 (FIB-4) — use qgcomp.boot for robust CI extraction
cat("\n[qgcomp] Stack 1 PFAS → fib4_log (qgcomp.boot B=200)...\n")
fit_qg <- tryCatch(
  qgcomp.boot(as.formula(paste(f_lhs, "~", f_rhs)),
              expnms = pfas_cols, data = df1, family = gaussian(),
              q = 4, B = 200, seed = 20260519, weights = wt_pooled),
  error = function(e) { cat("qgcomp.boot err:", conditionMessage(e), "\n"); NULL }
)
# fallback to noboot with safe extraction
if (is.null(fit_qg)) {
  cat("  fallback to qgcomp.noboot ...\n")
  fit_qg <- tryCatch(
    qgcomp.noboot(as.formula(paste(f_lhs, "~", f_rhs)),
                  expnms = pfas_cols, data = df1, family = gaussian(), q = 4, weights = wt_pooled),
    error = function(e) { cat("qgcomp.noboot err:", conditionMessage(e), "\n"); NULL }
  )
}
if (!is.null(fit_qg)) {
  # Robust psi/se extraction: try $psi/$var.psi (boot), then summary, then var.coef matrix
  psi <- NA_real_; psi_se <- NA_real_
  if (!is.null(fit_qg$psi) && !is.null(fit_qg$var.psi)) {
    psi <- as.numeric(fit_qg$psi)[1]
    vp <- as.numeric(fit_qg$var.psi)[1]
    psi_se <- if (!is.na(vp) && vp > 0) sqrt(vp) else NA_real_
  }
  if (is.na(psi)) {
    s_qg <- tryCatch(summary(fit_qg)$coefficients, error = function(e) NULL)
    if (!is.null(s_qg) && "psi1" %in% rownames(s_qg)) {
      psi <- s_qg["psi1","Estimate"]; psi_se <- s_qg["psi1","Std. Error"]
    } else if (!is.null(s_qg) && nrow(s_qg) >= 2) {
      psi <- s_qg[2,1]; psi_se <- s_qg[2,2]
    }
  }
  if (is.na(psi) && !is.null(fit_qg$coef) && !is.null(fit_qg$var.coef)) {
    vc <- fit_qg$var.coef
    psi <- as.numeric(fit_qg$coef)[2]
    psi_se <- if (is.matrix(vc) && nrow(vc) >= 2) sqrt(vc[2,2]) else sqrt(as.numeric(vc)[2])
  }
  if (!is.na(psi) && !is.na(psi_se)) {
    qg_result <- data.frame(
      method = "qgcomp", outcome = "fib4_log", n = nrow(df1),
      psi = psi, psi_low = psi - 1.96 * psi_se, psi_high = psi + 1.96 * psi_se,
      p = 2 * pnorm(-abs(psi / psi_se))
    )
    cat(sprintf("  qgcomp psi = %.4f (95%% CI: %.4f, %.4f), P = %.2e\n",
                psi, psi - 1.96 * psi_se, psi + 1.96 * psi_se,
                2 * pnorm(-abs(psi / psi_se))))
    saveRDS(fit_qg, "output/tables/qgcomp_stack1_fib4.rds")
  } else {
    cat("  WARN: qgcomp psi/se extraction failed; skipping.\n")
    qg_result <- data.frame()
  }
} else {
  qg_result <- data.frame()
}

# WQS Stack 1 — use only numeric/binary covariates to avoid 1-level factor issues
cat("\n[WQS] Stack 1 PFAS → fib4_log (numeric cov only)...\n")
fit_wqs <- tryCatch(
  gwqs(as.formula(paste(f_lhs, "~ wqs +", paste(cov_use, collapse = " + "))),
       mix_name = pfas_cols, data = df1, q = 4,
       validation = 0.3, family = "gaussian", b = 200,
       seed = 20260519, plots = FALSE),
  error = function(e) { cat("WQS err:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(fit_wqs)) {
  s <- summary(fit_wqs)$coefficients
  wqs_b <- as.numeric(s["wqs", 1]); wqs_se <- as.numeric(s["wqs", 2])
  # Robust n extraction
  n_wqs <- nrow(df1)
  pcol <- if ("Pr(>|t|)" %in% colnames(s)) "Pr(>|t|)" else colnames(s)[ncol(s)]
  wqs_result <- data.frame(
    method = "WQS", outcome = "fib4_log", n = n_wqs,
    psi = wqs_b, psi_low = wqs_b - 1.96 * wqs_se, psi_high = wqs_b + 1.96 * wqs_se,
    p = as.numeric(s["wqs", pcol]),
    stringsAsFactors = FALSE
  )
  cat(sprintf("  WQS beta = %.4f (95%% CI: %.4f, %.4f), P = %.2e\n",
              wqs_b, wqs_b - 1.96 * wqs_se, wqs_b + 1.96 * wqs_se,
              wqs_result$p))
  saveRDS(fit_wqs, "output/tables/wqs_stack1_fib4.rds")

  # WQS weights (defensive)
  wqs_weights <- tryCatch(fit_wqs$final_weights, error = function(e) NULL)
  if (!is.null(wqs_weights)) {
    write.csv(wqs_weights, "output/tables/wqs_weights_pfas.csv", row.names = FALSE)
  }
} else {
  wqs_result <- data.frame()
}

# Mixture triangulation summary
if (nrow(qg_result) + nrow(wqs_result) > 0) {
  triangulation <- rbind(qg_result, wqs_result)
  write.csv(triangulation, "output/tables/mixture_triangulation_stack1.csv", row.names = FALSE)
  cat("\nMixture triangulation (qgcomp + WQS) saved.\n")
}

cat("\n========================================\n")