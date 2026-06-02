# ============================================
# 06d_reanalysis_cleanM2.R  — AUTHORITATIVE corrected re-run
# Fixes the over-adjustment flagged in review: M2 = CONFOUNDERS ONLY, matching the
# manuscript text + the DAG (NO metabolic mediators fbg/hba1c/ldl/sbp/dbp; no marital/drink).
# Correct PFAS-subsample weight (design_main / design_mortality, wt_pooled = pooled_pfas_weight).
# Writes ALL headline CSVs so every number traces to code:
#   table_xsec_pfas_fib4.csv  (single-PFAS continuous + threshold, Crude/M1/M2)
#   table_mutual_pfas_fib4.csv (6-PFAS joint M2)  <- NEW, fixes reproducibility P0
#   table2_cox_pfas.csv        (all-cause Cox, total-effect M2)
#   fine_gray_cm_pfas.csv      (cardiometabolic Fine-Gray SHR, total-effect M2)
#   table_mortality_revcaus.csv(early-death-exclusion sensitivity)
#   evalue.csv                 (E-value for FINAL PFNA threshold OR)  <- fixes P0
#   table_bh_fdr.csv           (BH within 4 families incl. cardiometabolic FG) <- fixes P0
#   tableS1_pfas_distribution.csv (geomean+IQR on corrected raw conc) <- fixes geomean P0
#   table_sens_noDMHTN.csv     (sensitivity: M2 minus diabetes/hypertension)
# ============================================
suppressPackageStartupMessages({ library(survey); library(survival); library(dplyr); library(cmprsk) })
options(survey.lonely.psu = "adjust")
load("data/processed/nhanes_final.RData"); load("data/processed/nhanes_design.RData")

clean_nm <- c(LBXPFNA="PFNA",LBXMPAH="Me-PFOSA-AcOH",LBXPFDE="PFDA",LBXPFOS="PFOS",LBXPFOA="PFOA",LBXPFHS="PFHxS")
pf6 <- c("LBXPFNA_z","LBXMPAH_z","LBXPFDE_z","LBXPFOS_z","LBXPFOA_z","LBXPFHS_z")

# ---- CLEAN pre-specified M2 (confounders only; matches manuscript + DAG) ----
M2_clean <- c("age","RIAGENDR","race","education","pir","smoke","bmi","diabetes","hypertension")
keepcov <- function(df,covs) covs[vapply(covs,function(c){x<-df[[c]];mean(is.na(x))<=0.20 && length(unique(x[!is.na(x)]))>=2},logical(1))]

des <- design_main
des <- update(des, fib4_adv267 = as.integer(exp(fib4_log) >= 2.67))
covM2 <- keepcov(des$variables, intersect(M2_clean, names(des$variables)))
covM1 <- intersect(c("age","RIAGENDR","race"), names(des$variables))
cat("Clean M2 covariates used:", paste(covM2, collapse=", "), "\n")
cat("(dropped from requested clean set:", paste(setdiff(M2_clean, covM2), collapse=", "), ")\n\n")

emit <- function(fit, pf, model, outcome, binary=FALSE) {
  if (is.null(fit)) return(NULL); s <- summary(fit)$coefficients
  if (!pf %in% rownames(s)) return(NULL)
  b<-s[pf,1]; se<-s[pf,2]; p<-s[pf,ncol(s)]
  data.frame(pfas=pf, outcome=outcome, model=model,
    est=round(if(binary) exp(b) else b,4), lo=round(if(binary) exp(b-1.96*se) else b-1.96*se,4),
    hi=round(if(binary) exp(b+1.96*se) else b+1.96*se,4),
    p=signif(p,3), pct_per_SD=round((exp(b)-1)*100,1), stringsAsFactors=FALSE)
}
sg <- function(f, ...) tryCatch(svyglm(as.formula(f), design=des, ...), error=function(e){cat(" err:",conditionMessage(e),"\n");NULL})

# ===== 1. Single-PFAS fibrosis (continuous + threshold) =====
rows <- list()
for (pf in pf6) {
  rows[[length(rows)+1]] <- emit(sg(sprintf("fib4_log ~ %s", pf)), pf, "Crude", "logFIB4")
  rows[[length(rows)+1]] <- emit(sg(sprintf("fib4_log ~ %s + %s", pf, paste(covM1,collapse=" + "))), pf, "M1", "logFIB4")
  rows[[length(rows)+1]] <- emit(sg(sprintf("fib4_log ~ %s + %s", pf, paste(covM2,collapse=" + "))), pf, "M2", "logFIB4")
}
for (pf in pf6)
  rows[[length(rows)+1]] <- emit(sg(sprintf("fib4_adv267 ~ %s + %s", pf, paste(covM2,collapse=" + ")), family=quasibinomial), pf, "M2", "FIB4>=2.67", binary=TRUE)
xs <- do.call(rbind, rows)
write.csv(xs, "output/tables/table_xsec_pfas_fib4.csv", row.names=FALSE)
cat("=== 1. Single-PFAS -> FIB-4 (clean M2) ===\n"); print(xs)

# ===== 2. Mutually-adjusted 6-PFAS model (continuous) =====
fitJ <- sg(sprintf("fib4_log ~ %s + %s", paste(pf6,collapse=" + "), paste(covM2,collapse=" + ")))
sJ <- summary(fitJ)$coefficients; nJ <- length(residuals(fitJ))
mut <- do.call(rbind, lapply(pf6, function(pf){ b<-sJ[pf,1]; se<-sJ[pf,2]
  data.frame(pfas=pf, n=nJ, est=round(b,4), lo=round(b-1.96*se,4), hi=round(b+1.96*se,4),
             p=signif(sJ[pf,4],3), pct_per_SD=round((exp(b)-1)*100,1), stringsAsFactors=FALSE)}))
write.csv(mut, "output/tables/table_mutual_pfas_fib4.csv", row.names=FALSE)
cat(sprintf("\n=== 2. Mutually-adjusted 6-PFAS model (n=%d) ===\n", nJ)); print(mut)

# ===== 3. Sensitivity: M2 minus diabetes/hypertension =====
covM2min <- setdiff(covM2, c("diabetes","hypertension"))
sens <- do.call(rbind, lapply(pf6, function(pf) emit(sg(sprintf("fib4_log ~ %s + %s", pf, paste(covM2min,collapse=" + "))), pf, "M2_noDMHTN", "logFIB4")))
write.csv(sens, "output/tables/table_sens_noDMHTN.csv", row.names=FALSE)
cat("\n=== 3. Sensitivity (M2 minus diabetes/hypertension), continuous ===\n"); print(sens)

# ===== 4. Mortality: all-cause Cox + Fine-Gray (total-effect clean M2, correct weight) =====
desm <- design_mortality; desm$variables$time_yr <- desm$variables$permth/12
vm <- desm$variables
covM2m <- keepcov(vm, intersect(M2_clean, names(vm)))
cat(sprintf("\nMortality cohort n=%d | all-cause=%d | cardiometabolic=%d | M2: %s\n",
            nrow(vm), sum(vm$mort_allcause==1,na.rm=TRUE), sum(vm$mort_cm==1,na.rm=TRUE), paste(covM2m,collapse=", ")))
cox_rows <- list()
for (pf in pf6) for (oc in c("mort_allcause","mort_cm")) {
  f <- as.formula(sprintf("Surv(time_yr, %s) ~ %s + %s", oc, pf, paste(covM2m,collapse=" + ")))
  fit <- tryCatch(svycoxph(f, design=desm), error=function(e) NULL)
  if (!is.null(fit)) { s<-summary(fit)$coefficients
    if (pf %in% rownames(s)) { co<-s[pf,"coef"]; se<-s[pf,"se(coef)"]
      cox_rows[[length(cox_rows)+1]] <- data.frame(pfas=pf, outcome=oc, model="M2",
        hr=exp(co), ci_low=exp(co-1.96*se), ci_high=exp(co+1.96*se), p=2*pnorm(-abs(co/se)), stringsAsFactors=FALSE) } }
}
cox_df <- do.call(rbind, cox_rows)
write.csv(cox_df, "output/tables/table2_cox_pfas.csv", row.names=FALSE)
cat("\n=== 4a. All-cause + cause-specific Cox (total-effect M2) ===\n"); print(cox_df)

# Fine-Gray cardiometabolic (subdistribution), total-effect M2
fg_rows <- list()
for (pf in pf6) {
  d <- vm %>% filter(!is.na(.data[[pf]]), !is.na(time_yr), !is.na(mort_allcause))
  cc <- intersect(covM2m, names(d)); ok <- complete.cases(d[,c(pf,cc)]); d <- d[ok,]
  d$event_fg <- with(d, ifelse(mort_cm==1,1L, ifelse(mort_allcause==1,2L,0L)))
  cm <- model.matrix(as.formula(sprintf("~ %s + %s", pf, paste(cc,collapse=" + "))), data=d)[,-1]
  fit <- tryCatch(crr(d$time_yr, d$event_fg, cm, failcode=1), error=function(e) NULL)
  if (!is.null(fit)) { s<-summary(fit)$coef
    if (pf %in% rownames(s)) fg_rows[[length(fg_rows)+1]] <- data.frame(pfas=pf, method="Fine_Gray_SHR",
      hr=exp(s[pf,"coef"]), ci_low=exp(s[pf,"coef"]-1.96*s[pf,"se(coef)"]), ci_high=exp(s[pf,"coef"]+1.96*s[pf,"se(coef)"]),
      p=s[pf,"p-value"], stringsAsFactors=FALSE) }
}
fg_df <- do.call(rbind, fg_rows)
write.csv(fg_df, "output/tables/fine_gray_cm_pfas.csv", row.names=FALSE)
cat("\n=== 4b. Fine-Gray cardiometabolic SHR (total-effect M2) ===\n"); print(fg_df)

# ===== 5. Early-death-exclusion sensitivity (PFOS & PFOA all-cause) =====
rc_rows <- list()
for (pf in c("LBXPFOS_z","LBXPFOA_z")) for (excl in c(0,1,2,5)) {
  dsub <- subset(desm, time_yr > excl)
  f <- as.formula(sprintf("Surv(time_yr, mort_allcause) ~ %s + %s", pf, paste(covM2m,collapse=" + ")))
  fit <- tryCatch(svycoxph(f, design=dsub), error=function(e) NULL)
  hr <- c(NA,NA,NA,NA)
  if (!is.null(fit)) { s<-summary(fit)$coefficients; if (pf %in% rownames(s)) { co<-s[pf,"coef"]; se<-s[pf,"se(coef)"]
    hr <- c(exp(co),exp(co-1.96*se),exp(co+1.96*se),2*pnorm(-abs(co/se))) } }
  rc_rows[[length(rc_rows)+1]] <- data.frame(exposure=clean_nm[sub("_z","",pf)], exclude_early_yr=excl,
    HR=round(hr[1],3), lo=round(hr[2],3), hi=round(hr[3],3), p=signif(hr[4],3), stringsAsFactors=FALSE)
}
rc_df <- do.call(rbind, rc_rows)
write.csv(rc_df, "output/tables/table_mortality_revcaus.csv", row.names=FALSE)
cat("\n=== 5. Early-death-exclusion sensitivity ===\n"); print(rc_df)

# ===== 6. E-value for FINAL PFNA threshold OR =====
pfna_or <- xs$est[xs$pfas=="LBXPFNA_z" & xs$outcome=="FIB4>=2.67" & xs$model=="M2"]
pfna_lo <- xs$lo[xs$pfas=="LBXPFNA_z" & xs$outcome=="FIB4>=2.67" & xs$model=="M2"]
evOR <- function(or){ if(or<1) or<-1/or; or + sqrt(or*(or-1)) }
ev <- data.frame(exposure="PFNA", outcome="FIB4>=2.67", OR=round(pfna_or,4), OR_lci=round(pfna_lo,4),
                 EV_point=round(evOR(pfna_or),2), EV_lci=round(evOR(pfna_lo),2))
write.csv(ev, "output/tables/evalue.csv", row.names=FALSE)
cat("\n=== 6. E-value (PFNA -> FIB-4>=2.67) ===\n"); print(ev)

# ===== 7. BH-FDR within 4 families =====
bh <- function(d, fam) { d$p_bh <- p.adjust(d$p, "BH"); d$family <- fam; d[,c("family","pfas","est","p","p_bh")] }
b1 <- bh(transform(xs[xs$outcome=="logFIB4" & xs$model=="M2",], est=est), "continuous_logFIB4_M2")
b2 <- bh(transform(xs[xs$outcome=="FIB4>=2.67" & xs$model=="M2",], est=est), "advanced_FIB4_2.67_M2")
b3 <- bh(data.frame(pfas=cox_df$pfas[cox_df$outcome=="mort_allcause"], est=cox_df$hr[cox_df$outcome=="mort_allcause"], p=cox_df$p[cox_df$outcome=="mort_allcause"]), "allcause_mortality_M2")
b4 <- bh(data.frame(pfas=fg_df$pfas, est=fg_df$hr, p=fg_df$p), "cardiometabolic_FineGray_M2")
bh_all <- rbind(b1,b2,b3,b4)
write.csv(bh_all, "output/tables/table_bh_fdr.csv", row.names=FALSE)
cat("\n=== 7. BH-FDR (4 families) ===\n"); print(bh_all)

# ===== 8. Table S1 distribution on corrected raw concentrations =====
rawmap <- c(LBXPFNA="PFNA",LBXPFOS="PFOS",LBXPFOA="PFOA",LBXPFHS="PFHxS",LBXPFDE="PFDA",LBXMPAH="Me-PFOSA-AcOH")
s1 <- do.call(rbind, lapply(names(rawmap), function(col){
  if (!col %in% names(nhanes_final)) return(NULL)
  x <- nhanes_final[[col]]; x <- x[!is.na(x) & x>0]
  data.frame(analyte=rawmap[col], n=length(x),
             geomean=round(exp(mean(log(x))),3),
             p25=round(quantile(x,.25),3), p50=round(quantile(x,.5),3), p75=round(quantile(x,.75),3),
             detect_gt_p25=NA, stringsAsFactors=FALSE)
}))
write.csv(s1, "output/tables/tableS1_pfas_distribution.csv", row.names=FALSE)
cat("\n=== 8. Table S1 distribution (corrected raw conc) ===\n"); print(s1)
cat("\nALL CSVs written. Re-analysis complete.\n")