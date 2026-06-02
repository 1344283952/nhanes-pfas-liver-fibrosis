# — age-free triangulation for the null (addresses FIB-4 age-in-numerator concern).
# APRI = (AST/ULN)/platelet x 100 has NO age term, so if PFAS are null vs APRI (with age still
# adjusted as a confounder), the absence of an independent association is not a FIB-4-algebra artifact.
suppressPackageStartupMessages({ library(survey); library(dplyr) })
load("data/processed/nhanes_final.RData"); load("data/processed/nhanes_design.RData")
des <- design_main
stopifnot("apri" %in% names(des$variables))
des <- update(des, log_apri = log(pmax(apri, 0.01)))
M2 <- intersect(c("age","RIAGENDR","race","education","pir","smoke","bmi","diabetes","hypertension"), names(des$variables))
pf6 <- c("LBXPFNA_z","LBXPFOS_z","LBXPFOA_z","LBXPFHS_z","LBXPFDE_z","LBXMPAH_z")
nm  <- c(LBXPFNA_z="PFNA",LBXPFOS_z="PFOS",LBXPFOA_z="PFOA",LBXPFHS_z="PFHxS",LBXPFDE_z="PFDA",LBXMPAH_z="Me-PFOSA-AcOH")
fit1 <- function(pf, rhs){
  f <- as.formula(paste("log_apri ~", paste(c(pf,rhs), collapse=" + ")))
  s <- summary(svyglm(f, design=des))$coef[pf,]
  c(pct=unname((exp(s["Estimate"])-1)*100), p=unname(s["Pr(>|t|)"]))
}
out <- data.frame(pfas=nm[pf6],
  crude_pct=sapply(pf6,function(p) round(fit1(p,character(0))["pct"],1)),
  crude_p =sapply(pf6,function(p) signif(fit1(p,character(0))["p"],3)),
  m2_pct  =sapply(pf6,function(p) round(fit1(p,M2)["pct"],1)),
  m2_p    =sapply(pf6,function(p) signif(fit1(p,M2)["p"],3)), row.names=NULL)
write.csv(out, "output/tables/table_agefree_apri.csv", row.names=FALSE)
cat("=== PFAS vs log(APRI) [age-free fibrosis index] ===\n"); print(out)
cat(sprintf("\nAPRI available n = %d\n", sum(!is.na(des$variables$apri))))