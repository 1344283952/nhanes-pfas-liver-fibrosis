# ============================================
# scripts/04_survey_design.R
# 复杂抽样设计 (NHANES weighted, pooled C-I + P_ per NCHS Series 2 No. 190)
# 输入: data/processed/nhanes_final.RData
# 输出: data/processed/nhanes_design.RData
#
# : 3 个 design 对应 3 个 cohort stack:
#   design_main      → Stack 1 主分析 (N≈12,533)
#   design_fibroscan → Stack 2 FibroScan 子 (N≈2,391, P_ only)
#   design_metals    → Stack 4 PFAS+Metals (N≈2,500+)
# ============================================

library(survey); library(dplyr)

cat("========================================\n")
cat("复杂抽样设计 (3 stack)\n")
cat("========================================\n\n")

load("data/processed/nhanes_final.RData")
options(survey.lonely.psu = "adjust")

# ----------------------------------------------------------
# Design 1: Stack 1 主分析 (PFAS D-I + P_, J⊂P_ dedupe)
# ----------------------------------------------------------
design_main <- svydesign(
  ids = ~SDMVPSU, strata = ~SDMVSTRA,
  weights = ~wt_pooled, data = nhanes_final,
  nest = TRUE
)
cat(sprintf("design_main (Stack 1 主): n=%d, weights=wt_pooled (NCHS Series 2 No. 190)\n",
            nrow(nhanes_final)))

# Subset designs (用 subset 而不是 svydesign with subset df, 保留抽样信息)
design_mortality <- subset(design_main, !is.na(ELIGSTAT) & ELIGSTAT == 1)
cat(sprintf("design_mortality (Stack 7 cardiometabolic): n=%d\n",
            nrow(design_mortality$variables)))

# ----------------------------------------------------------
# Design 2: Stack 2 FibroScan (P_ only)
# ----------------------------------------------------------
# P_ FibroScan design: PFAS-subsample serum weight WTSBAPRP (consistent with all PFAS estimates), NOT MEC
if (nrow(df_fibroscan) > 0) {
  design_fibroscan <- svydesign(
    ids = ~SDMVPSU, strata = ~SDMVSTRA,
    weights = ~WTSBAPRP, data = df_fibroscan,
    nest = TRUE
  )
  cat(sprintf("design_fibroscan (Stack 2 FibroScan P_): n=%d, weights=WTSBAPRP\n",
              nrow(df_fibroscan)))
} else {
  design_fibroscan <- NULL
  cat("design_fibroscan: NULL (df_fibroscan 为空)\n")
}

# ----------------------------------------------------------
# Design 3: Stack 4 PFAS+Metals (复用 wt_pooled)
# ----------------------------------------------------------
if (nrow(df_metals) > 0) {
  design_metals <- svydesign(
    ids = ~SDMVPSU, strata = ~SDMVSTRA,
    weights = ~wt_pooled, data = df_metals,
    nest = TRUE
  )
  cat(sprintf("design_metals (Stack 4 PFAS+Metals): n=%d, weights=wt_pooled\n",
              nrow(df_metals)))
} else {
  design_metals <- NULL
}

# Save
save(design_main, design_mortality, design_fibroscan, design_metals,
     file = "data/processed/nhanes_design.RData")
cat("\n已保存 data/processed/nhanes_design.RData (4 design)\n")
cat("========================================\n")