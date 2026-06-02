# ============================================
# scripts/_pfas_weight.R  (project-local; templates/ stay read-only)
# PFAS one-third-subsample POOLED weight, per NCHS Series 2 No. 190.
#
# WHY: serum PFAS were measured in a ONE-THIRD subsample of persons >=12y; the CDC
# documentation states "special sample weights are required". The full MEC weight
# (WTMEC2YR/WTMECPRP) over-counts the 2/3 not in the PFAS subsample and biases every
# weighted estimate. Correct base weight = the PFAS subsample serum weight.
#
# The PFAS subsample-weight variable name DRIFTS across cycles (verified empirically
# from the merged data: each column fully covers its cycle's PFAS-measured rows, and
# confirmed against CDC PFAS_H.htm = WTSB2YR / P_PFAS.htm = WTSBAPRP):
#   WTSA2YR  -> 2005-06 (D), 2011-12 (G)
#   WTSC2YR  -> 2007-08 (E), 2009-10 (F)
#   WTSB2YR  -> 2013-14 (H), 2015-16 (I)
#   WTSBAPRP -> 2017-March 2020 pre-pandemic (P_)
#
# Pooling (NCHS Series 2 No.190 §4): weight x (cycle_years / total_years).
# ============================================

pooled_pfas_weight <- function(df, cycle_years) {
  total_years <- sum(cycle_years)
  yrs   <- cycle_years[as.character(df$cycle_tag)]
  is_pp <- grepl("PrePandemic|^P_", as.character(df$cycle_tag))
  # 2-yr cycles: coalesce the (mutually exclusive) per-cycle PFAS subsample weights
  base2 <- ifelse(!is.na(df$WTSA2YR) & df$WTSA2YR > 0, df$WTSA2YR,
           ifelse(!is.na(df$WTSC2YR) & df$WTSC2YR > 0, df$WTSC2YR, df$WTSB2YR))
  base  <- ifelse(is_pp, df$WTSBAPRP, base2)
  base * (yrs / total_years)
}