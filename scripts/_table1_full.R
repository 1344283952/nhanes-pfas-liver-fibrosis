# — full Table 1 by PFNA tertile (M2 covariates) -> table1_by_pfna.csv (traceable)
suppressPackageStartupMessages({ library(survey); library(dplyr) })
load("data/processed/nhanes_final.RData"); load("data/processed/nhanes_design.RData")
des <- update(design_main, adv = as.integer(exp(fib4_log) >= 2.67),
  pfna_t = cut(LBXPFNA_z, quantile(LBXPFNA_z, c(0,1/3,2/3,1), na.rm=TRUE), labels=c("T1","T2","T3"), include.lowest=TRUE))
v <- des$variables
rows <- list()
add <- function(lab, t1, t2, t3) rows[[length(rows)+1]] <<- data.frame(variable=lab, T1=t1, T2=t2, T3=t3, stringsAsFactors=FALSE)
numrow <- function(vn, lab){ m <- svyby(as.formula(paste0("~",vn)), ~pfna_t, des, svymean, na.rm=TRUE); add(lab, sprintf("%.1f",m[[vn]][1]), sprintf("%.1f",m[[vn]][2]), sprintf("%.1f",m[[vn]][3])) }
pctrow <- function(cond, lab){ des2 <- update(des, .x=as.integer(cond)); m <- svyby(~.x, ~pfna_t, des2, svymean, na.rm=TRUE); add(lab, sprintf("%.1f",m$.x[1]*100), sprintf("%.1f",m$.x[2]*100), sprintf("%.1f",m$.x[3]*100)) }

nt <- table(v$pfna_t); add("n (unweighted)", nt[1], nt[2], nt[3])
numrow("age", "Age, years (mean)")
pctrow(v$RIAGENDR==2, "Female, %")
if ("race" %in% names(v)) {
  pctrow(as.character(v$race)=="Non-Hispanic White", "Race: Non-Hispanic White, %")
  pctrow(as.character(v$race)=="Non-Hispanic Black", "Race: Non-Hispanic Black, %")
  pctrow(as.character(v$race)=="Mexican American", "Race: Mexican American, %")
  pctrow(as.character(v$race) %in% c("Other Hispanic","Other Race"), "Race: other, %")
}
if ("education" %in% names(v)) pctrow(as.character(v$education)=="College or above", "College or above, %")
numrow("pir", "Poverty-income ratio (mean)")
if ("smoke" %in% names(v)) pctrow(as.character(v$smoke)=="Ever", "Ever smoker, %")
numrow("bmi", "Body-mass index, kg/m^2 (mean)")
if ("diabetes" %in% names(v)) pctrow(as.character(v$diabetes) %in% c("Yes","1"), "Diabetes, %")
if ("hypertension" %in% names(v)) pctrow(as.character(v$hypertension) %in% c("Yes","1"), "Hypertension, %")
# family history of CVD removed: raw NHANES variable was uncleaned (codes 1/2/7/9) and not in the M2 model
if ("LBXSATSI" %in% names(v)) numrow("LBXSATSI", "ALT, U/L (mean)")
if ("LBXSASSI" %in% names(v)) numrow("LBXSASSI", "AST, U/L (mean)")
pctrow(v$adv==1, "FIB-4 >= 2.67, %")

t1 <- do.call(rbind, rows)
write.csv(t1, "output/tables/table1_by_pfna.csv", row.names=FALSE)
cat("Wrote expanded table1_by_pfna.csv:\n"); print(t1)
cat(sprintf("\nFIB-4>=2.67 events (unweighted): %d total (%s by tertile)\n",
    sum(v$adv,na.rm=TRUE), paste(tapply(v$adv,v$pfna_t,sum,na.rm=TRUE),collapse="/")))