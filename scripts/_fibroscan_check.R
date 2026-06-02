# — gold-standard age-free triangulation: PFAS vs FibroScan liver stiffness (LSM, kPa).
# LSM is a direct elastography fibrosis measure with NO age term, so it adjudicates the
# FIB-4 (age-embedded, null) vs APRI (age-free, positive) discrepancy.
suppressPackageStartupMessages({ library(survey); library(dplyr) })
load("data/processed/nhanes_final.RData"); load("data/processed/nhanes_design.RData")
options(survey.lonely.psu = "adjust")
# Use the PFAS-subsample serum weight (WTSBAPRP) — consistent with every other PFAS estimate —
# rather than the full MEC examination weight, on the FibroScan subset.
fv <- design_fibroscan$variables
des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTSBAPRP, data = fv, nest = TRUE)
v <- des$variables
cat(sprintf("FibroScan subset n=%d; has lsm_median=%s; PFNA non-NA=%d; CAP non-NA=%d\n",
            nrow(v), "lsm_median" %in% names(v), sum(!is.na(v$LBXPFNA_z)),
            if ("cap_median" %in% names(v)) sum(!is.na(v$cap_median)) else 0))
des <- update(des, log_lsm = log(pmax(lsm_median, 0.1)))
M2 <- intersect(c("age","RIAGENDR","race","education","pir","smoke","bmi","diabetes","hypertension"), names(v))
pf6 <- c("LBXPFNA_z","LBXPFOS_z","LBXPFOA_z","LBXPFHS_z","LBXPFDE_z","LBXMPAH_z")
nm  <- c(LBXPFNA_z="PFNA",LBXPFOS_z="PFOS",LBXPFOA_z="PFOA",LBXPFHS_z="PFHxS",LBXPFDE_z="PFDA",LBXMPAH_z="Me-PFOSA-AcOH")
fit1 <- function(pf, rhs){
  f <- as.formula(paste("log_lsm ~", paste(c(pf,rhs), collapse=" + ")))
  m <- tryCatch(svyglm(f, design=des), error=function(e) NULL)
  if (is.null(m) || !(pf %in% rownames(summary(m)$coef))) return(c(pct=NA,p=NA))
  s <- summary(m)$coef[pf,]; c(pct=unname((exp(s["Estimate"])-1)*100), p=unname(s["Pr(>|t|)"]))
}
out <- data.frame(pfas=nm[pf6],
  crude_pct=sapply(pf6,function(p) round(fit1(p,character(0))["pct"],1)),
  crude_p =sapply(pf6,function(p) signif(fit1(p,character(0))["p"],3)),
  m2_pct  =sapply(pf6,function(p) round(fit1(p,M2)["pct"],1)),
  m2_p    =sapply(pf6,function(p) signif(fit1(p,M2)["p"],3)), row.names=NULL)
write.csv(out, "output/tables/table_agefree_lsm.csv", row.names=FALSE)
cat("=== PFAS vs log(FibroScan LSM) [gold-standard age-free fibrosis, M2 incl. age] ===\n"); print(out)