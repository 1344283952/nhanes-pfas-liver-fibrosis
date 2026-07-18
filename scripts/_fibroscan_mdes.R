# — 弹性成像(LSM)子样本的最小可检测效应 MDES
# 目的: 审稿人必问 "n=1,929 的 null 是不是功效不足?" —— 用数字回答, 不用形容词。
# MDES = (z_{0.975} + z_{0.80}) * SE = 2.802 * SE, 即 80% 功效下能检出的最小效应。
# 与 _fibroscan_check.R 同口径(同设计/同权重/同 M2), 只额外取 SE。
suppressPackageStartupMessages({ library(survey); library(dplyr) })
load("data/processed/nhanes_final.RData"); load("data/processed/nhanes_design.RData")
options(survey.lonely.psu = "adjust")

fv  <- design_fibroscan$variables
des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTSBAPRP, data = fv, nest = TRUE)
des <- update(des, log_lsm = log(pmax(lsm_median, 0.1)))
v   <- des$variables

M2  <- intersect(c("age","RIAGENDR","race","education","pir","smoke","bmi","diabetes","hypertension"), names(v))
pf6 <- c("LBXPFNA_z","LBXPFOS_z","LBXPFOA_z","LBXPFHS_z","LBXPFDE_z","LBXMPAH_z")
nm  <- c(LBXPFNA_z="PFNA", LBXPFOS_z="PFOS", LBXPFOA_z="PFOA",
         LBXPFHS_z="PFHxS", LBXPFDE_z="PFDA", LBXMPAH_z="Me-PFOSA-AcOH")

K <- qnorm(0.975) + qnorm(0.80)   # = 2.8016
rows <- list()
for (pf in pf6) {
  f <- as.formula(paste("log_lsm ~", paste(c(pf, M2), collapse = " + ")))
  m <- tryCatch(svyglm(f, design = des), error = function(e) NULL)
  if (is.null(m)) next
  s  <- summary(m)$coef[pf, ]
  se <- unname(s["Std. Error"])
  rows[[length(rows)+1]] <- data.frame(
    pfas          = nm[[pf]],
    n_model       = nobs(m),
    obs_pct_perSD = round((exp(unname(s["Estimate"])) - 1) * 100, 2),
    p             = signif(unname(s["Pr(>|t|)"]), 3),
    se_log        = round(se, 5),
    mdes_pct_perSD= round((exp(K * se) - 1) * 100, 2),   # 80% 功效可检出的最小 % 变化
    stringsAsFactors = FALSE)
}
out <- do.call(rbind, rows)
write.csv(out, "output/tables/table_lsm_mdes.csv", row.names = FALSE)
cat("=== 弹性成像 LSM 子样本 MDES (80% 功效, 双侧 α=0.05) ===\n"); print(out)
cat(sprintf("\nMDES 范围: %.1f%% – %.1f%% / SD;  实测最大|效应| = %.1f%%\n",
            min(out$mdes_pct_perSD), max(out$mdes_pct_perSD), max(abs(out$obs_pct_perSD))))