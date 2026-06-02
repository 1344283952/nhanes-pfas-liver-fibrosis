# ============================================
# scripts/03_clean_data.R
# 输入: data/processed/nhanes_raw_merged.RData (nhanes_all 9 周期 + rx_all + mort_all)
# 输出: data/processed/nhanes_final.RData (nhanes_final + df_pp_subset + df_fibroscan)
#
# 23 加项关键 logic:
#   清洗加项: HEPB/HEPC 排除 + Cotinine + Fish freq + RHQ postmenopausal + HSCRP
#   PFAS 处理: LOD/√2 substitution + log2-transform + z-score + Σ-PFAS-C/-S/mole-ratio
#   5-score: FIB-4 + AGILE 3+ + AGILE 4 + FAST + APRI + Hepamet HFS
#   MetALD 三组 (Rinella 2023): 复用 templates/_shared/metald_classify.R
#   COV_M2_COX (templates/_shared/covariates_config.R) + 客观吸烟 (LBXCOT) + fish + drinking water
#
# J⊂P_ dedupe: 用 cycle_tag, 主分析排除 NHANES_2017_2018 (J 行)
# Cohort: Stack 1 主 N=12,533 / Stack 2 FibroScan P_ N=2,391 / Stack 4 PFAS+Metals N~2,545+
# ============================================

library(dplyr); library(tidyr); library(purrr)

source("templates/_shared/fib4.R")
source("templates/_shared/hepamet_hfs.R")
source("templates/_shared/metald_classify.R")
source("templates/_shared/pooled_weight.R")
source("templates/_shared/var_aliases.R")
source("templates/_shared/cycle_config.R")
source("templates/_shared/diabetes_define.R")
source("templates/_shared/covariates_config.R")

set.seed(20260519)  # 启动日 reproducibility seed

cat("========================================\n")
cat("清洗 + 23 加项 (Stack 1 主 N≈12,533)\n")
cat("========================================\n\n")

load("data/processed/nhanes_raw_merged.RData")
cat(sprintf("原始合并: %d 行 x %d 列\n\n", nrow(nhanes_all), ncol(nhanes_all)))

flow <- list()
log_flow <- function(label, n) {
  cat(sprintf("  [流程] %-60s n = %d\n", label, n))
  flow[[length(flow) + 1]] <<- data.frame(step = label, n = n, stringsAsFactors = FALSE)
}
log_flow("原始合并 (9 cycles, C-J + P_)", nrow(nhanes_all))

# Util
na_codes <- function(x, codes) ifelse(x %in% codes, NA, x)
zero_to_na <- function(x) ifelse(!is.na(x) & x == 0, NA, x)

# Step 0: 跨周期变量名 normalize (var_aliases.R)
nhanes_all <- apply_var_aliases(nhanes_all)
cat("apply_var_aliases: canonical 列已生成\n\n")

# FIX 2026-06-01 ( P0): templates/_shared/var_aliases.R transposed AST/ALT.
# CDC BIOPRO codebook (verified 2026-06-01): LBXSATSI = "Alanine Aminotransferase (ALT) (U/L)";
# LBXSASSI = "Aspartate Aminotransferase (AST) (U/L)". The alias map set ast_unl<-LBXSATSI (ALT)
# and alt_unl<-LBXSASSI (AST), swapping the two enzymes. Reassign from the correctly-identified
# raw columns so FIB-4 = (age*AST)/(PLT*sqrt(ALT)) per Sterling 2006. apply_var_aliases() keeps
# the raw LBX* columns, so LBXSATSI/LBXSASSI are still present here.
nhanes_all$ast_unl <- coalesce_cols(nhanes_all, c("LBXSASSI","LBXSAST"))  # AST (Aspartate)
nhanes_all$alt_unl <- coalesce_cols(nhanes_all, c("LBXSATSI","LBXSALT"))  # ALT (Alanine)
cat(sprintf("AST/ALT swap fix applied: mean AST=%.1f U/L, mean ALT=%.1f U/L (n AST=%d, ALT=%d)\n\n",
            mean(nhanes_all$ast_unl, na.rm=TRUE), mean(nhanes_all$alt_unl, na.rm=TRUE),
            sum(!is.na(nhanes_all$ast_unl)), sum(!is.na(nhanes_all$alt_unl))))

# ----------------------------------------------------------
# Step 1: 入选 — age ≥ 20 + 非妊娠
# ----------------------------------------------------------
df <- nhanes_all %>% filter(RIDAGEYR >= 20)
log_flow("Age >= 20", nrow(df))
df <- df %>% filter(is.na(RIDEXPRG) | RIDEXPRG != 1)
log_flow("+ not pregnant (excl confirmed RIDEXPRG==1 only)", nrow(df))

# ----------------------------------------------------------
# Step 2: PFAS 6 核心化合物 任一非缺 (Stack 1 主入选)
# 跨周期 PFAS 字段名: LBXPFOA / LBDPFOA (检测 LOD 标志); LBXPFOS / LBDPFOS; etc.
# 实际 NHANES 字段: LBXPFOA, LBXPFOS, LBXPFNA, LBXPFHS (PFHxS), LBXPFDE (PFDA), LBXMPAH (Me-PFOSA-AcOH)
# ----------------------------------------------------------
# FIX (2026-06-01): reconstruct TOTAL PFOS/PFOA from linear+branched isomers where the
# combined column is absent. 2005-2012 (PFC_*) report combined LBXPFOS/LBXPFOA; 2013-14
# (SSPFAS_H) report SSNPFOS+SSMPFOS / SSNPFOA+SSBPFOA; 2015-March 2020 (PFAS_I/J/P_) report
# LBXNFOS+LBXMFOS / LBXNFOA+LBXBFOA. Without this, PFOS/PFOA were NA for 49% of the cohort
# (all 2013-2020 cycles) and floor-imputed, making the "legacy null" a mechanical artifact.
.sum2 <- function(a, b) { x<-ifelse(is.na(a),0,a); y<-ifelse(is.na(b),0,b); o<-x+y; o[is.na(a)&is.na(b)]<-NA; o }
if (!"LBXPFOS" %in% names(df)) df$LBXPFOS <- NA_real_
if (!"LBXPFOA" %in% names(df)) df$LBXPFOA <- NA_real_
if (all(c("LBXNFOS","LBXMFOS") %in% names(df))) df$LBXPFOS <- ifelse(is.na(df$LBXPFOS), .sum2(df$LBXNFOS, df$LBXMFOS), df$LBXPFOS)
if (all(c("LBXNFOA","LBXBFOA") %in% names(df))) df$LBXPFOA <- ifelse(is.na(df$LBXPFOA), .sum2(df$LBXNFOA, df$LBXBFOA), df$LBXPFOA)
if (all(c("SSNPFOS","SSMPFOS") %in% names(df))) df$LBXPFOS <- ifelse(is.na(df$LBXPFOS), .sum2(df$SSNPFOS, df$SSMPFOS), df$LBXPFOS)
if (all(c("SSNPFOA","SSBPFOA") %in% names(df))) df$LBXPFOA <- ifelse(is.na(df$LBXPFOA), .sum2(df$SSNPFOA, df$SSBPFOA), df$LBXPFOA)
cat(sprintf("Total PFOS non-NA: %d | PFOA non-NA: %d (post isomer reconstruction)\n",
            sum(!is.na(df$LBXPFOS)), sum(!is.na(df$LBXPFOA))))

pfas_cols_core <- c("LBXPFOA","LBXPFOS","LBXPFNA","LBXPFHS","LBXPFDE","LBXMPAH")
# 兼容部分周期字段名漂移: 若缺则用 alternate name
pfas_alt_map <- list(
  "LBXPFHS" = c("LBXPFHS","LBXPFHxS"),
  "LBXPFDE" = c("LBXPFDE","LBXPFDA"),
  "LBXMPAH" = c("LBXMPAH","LBXMEPA")
)
for (canon in names(pfas_alt_map)) {
  if (!canon %in% names(df)) {
    for (alt in pfas_alt_map[[canon]][-1]) {
      if (alt %in% names(df)) {
        df[[canon]] <- df[[alt]]
        break
      }
    }
  }
}
# 实际可用 PFAS 列 (avoid all-NA cols)
pfas_cols_avail <- intersect(pfas_cols_core, names(df))
pfas_cols_avail <- pfas_cols_avail[sapply(pfas_cols_avail, function(c) sum(!is.na(df[[c]])) > 0)]
cat(sprintf("PFAS 核心化合物可用 (有数据): %s\n", paste(pfas_cols_avail, collapse = ", ")))

df <- df %>% filter(rowSums(!is.na(across(all_of(pfas_cols_avail)))) > 0)
log_flow(sprintf("+ PFAS 核心 (%d 化合物) 任一非缺", length(pfas_cols_avail)), nrow(df))

# ----------------------------------------------------------
# Step 3: PFAS LOD/√2 substitution + log2-transform + z-score
# LOD 标志列: LBDPFOA (1 = 检出, 0 = <LOD)
# 简化处理: NA → 视为 <LOD, 用 LOD/√2 = (各化合物 5th percentile) / sqrt(2)
# ----------------------------------------------------------
for (col in pfas_cols_avail) {
  vals <- df[[col]]
  # FIX: PFOS/PFOA in 2013-14 were measured ONLY in a separate surplus subsample
  # (SSPFAS_H, weight WTSSBH2Y) that is not poolable with the main PFAS subsample -> treat as
  # NOT-MEASURED and exclude (do NOT floor-impute). PFOS/PFOA therefore use six cycles
  # (2005-2012 + 2015-March 2020); PFNA and the other analytes use all seven.
  notmeas <- if (col %in% c("LBXPFOS","LBXPFOA")) df$cycle_tag == "NHANES_2013_2014" else rep(FALSE, length(vals))
  vals[notmeas] <- NA
  if (sum(!is.na(vals)) < 50) next
  lod_est <- min(vals[vals > 0], na.rm = TRUE)
  imp_val <- lod_est / sqrt(2)
  imp <- ifelse(is.na(vals) | vals <= 0, imp_val, vals)
  imp[notmeas] <- NA   # exclude not-measured rows (no floor imputation)
  df[[paste0(col, "_imp")]]  <- imp
  df[[paste0(col, "_log2")]] <- log2(imp)
  df[[paste0(col, "_z")]]    <- as.numeric(scale(log2(imp)))
}

# Σ-PFAS-C (carboxylic): PFOA + PFNA + PFDA + PFUnDA + PFDoDA (用可得的)
sum_C_cols <- intersect(c("LBXPFOA_imp","LBXPFNA_imp","LBXPFDE_imp"), names(df))
if (length(sum_C_cols) >= 2) {
  df$sum_pfas_c <- rowSums(df[, sum_C_cols], na.rm = TRUE)
  df$sum_pfas_c_log2 <- log2(pmax(df$sum_pfas_c, 0.01))
  df$sum_pfas_c_z <- as.numeric(scale(df$sum_pfas_c_log2))
}
# Σ-PFAS-S (sulfonic): PFOS + PFHxS
sum_S_cols <- intersect(c("LBXPFOS_imp","LBXPFHS_imp"), names(df))
if (length(sum_S_cols) >= 2) {
  df$sum_pfas_s <- rowSums(df[, sum_S_cols], na.rm = TRUE)
  df$sum_pfas_s_log2 <- log2(pmax(df$sum_pfas_s, 0.01))
  df$sum_pfas_s_z <- as.numeric(scale(df$sum_pfas_s_log2))
}
# PFOA/PFOS mole ratio (用 imp 值)
if ("LBXPFOA_imp" %in% names(df) && "LBXPFOS_imp" %in% names(df)) {
  # PFOA MW=414, PFOS MW=500
  df$pfoa_pfos_mole_ratio <- (df$LBXPFOA_imp / 414) / (df$LBXPFOS_imp / 500)
  df$pfoa_pfos_mole_ratio[!is.finite(df$pfoa_pfos_mole_ratio)] <- NA
}
cat(sprintf("PFAS 衍生暴露: sum_pfas_c (n=%d), sum_pfas_s (n=%d), pfoa_pfos_mole_ratio (n=%d)\n",
            sum(!is.na(df$sum_pfas_c %||% rep(NA, nrow(df)))),
            sum(!is.na(df$sum_pfas_s %||% rep(NA, nrow(df)))),
            sum(!is.na(df$pfoa_pfos_mole_ratio %||% rep(NA, nrow(df))))))

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ----------------------------------------------------------
# Step 4: HEPB/HEPC 排除 (Top #1, 必加不加 desk reject)
# ----------------------------------------------------------
df$hep_excl <- (df$hbsag %in% 1) | (df$hcv %in% 1)
df$hep_excl[is.na(df$hep_excl)] <- FALSE
n_hep_pos <- sum(df$hep_excl, na.rm = TRUE)
df <- df %>% filter(!hep_excl)
log_flow(sprintf("+ 排除 HEPB/HEPC 阳性 (n=%d)", n_hep_pos), nrow(df))

# ----------------------------------------------------------
# Step 5: 排除 self-reported 已确诊肝病 (MCQ160L == 1)
# ----------------------------------------------------------
df$MCQ160L <- na_codes(df$MCQ160L, c(7, 9))
n_liver_pre <- sum(df$MCQ160L == 1, na.rm = TRUE)
df <- df %>% filter(!(MCQ160L %in% 1))
log_flow(sprintf("+ 排除 self-reported 已诊断肝病 (n=%d)", n_liver_pre), nrow(df))

# ----------------------------------------------------------
# Step 6: FIB-4 + 衍生 (复用 _shared/fib4.R)
# ----------------------------------------------------------
df <- df %>% filter(!is.na(ast_unl), !is.na(alt_unl), !is.na(LBXPLTSI))
log_flow("+ FIB-4 ready (ALT/AST/PLT 非缺)", nrow(df))

df$fib4 <- calc_fib4(df$RIDAGEYR, df$ast_unl, df$alt_unl, df$LBXPLTSI)
df$fib4_log <- log(pmax(df$fib4, 0.01))
df$fib4_advanced <- as.integer(df$fib4 >= 2.67)
df$fib4_primary <- fib4_advanced_age_adj(df$fib4, df$RIDAGEYR)
df$fib4_high <- fib4_high_risk(df$fib4)

# APRI = (AST / ULN) / Platelet × 100; ULN = 40 U/L (women) 50 (men) 简化 40
df$apri <- ((df$ast_unl / 40) / df$LBXPLTSI) * 100
df$apri_advanced <- as.integer(df$apri >= 1.0)
df$apri_significant <- as.integer(df$apri >= 0.5)

# ----------------------------------------------------------
# Step 7: 血压均值 + HOMA-IR + HFS (Hepamet)
# ----------------------------------------------------------
df <- df %>% mutate(
  BPXSY1 = zero_to_na(BPXSY1), BPXSY2 = zero_to_na(BPXSY2),
  BPXSY3 = zero_to_na(BPXSY3), BPXSY4 = zero_to_na(BPXSY4),
  BPXDI1 = zero_to_na(BPXDI1), BPXDI2 = zero_to_na(BPXDI2),
  BPXDI3 = zero_to_na(BPXDI3), BPXDI4 = zero_to_na(BPXDI4)
) %>% rowwise() %>% mutate(
  sbp = mean(c(BPXSY2, BPXSY3, BPXSY4), na.rm = TRUE),
  dbp = mean(c(BPXDI2, BPXDI3, BPXDI4), na.rm = TRUE)
) %>% ungroup() %>% mutate(
  sbp = ifelse(is.nan(sbp), NA, sbp),
  dbp = ifelse(is.nan(dbp), NA, dbp)
)

df <- df %>% mutate(
  homa_ir = ifelse(!is.na(LBXIN) & !is.na(LBXGLU), (LBXIN * LBXGLU) / 405, NA),
  sex_male = as.integer(RIAGENDR == 1),
  albumin_gdl = alb_gl / 10
)

# Diabetes (ADA 2025 7-criteria via _shared/diabetes_define.R)
if (exists("rx_all") && nrow(rx_all) > 0 && "RXDDRUG" %in% names(rx_all)) {
  rx_dm <- detect_antidiabetic_rx(rx_all)
  df <- df %>% left_join(rx_dm, by = "SEQN") %>%
    mutate(rx_diabetes_yes = ifelse(is.na(rx_diabetes_yes), FALSE, rx_diabetes_yes))
}
df$diabetes_raw <- diabetes_comprehensive(df)
df$diabetes <- as.integer(ifelse(is.na(df$diabetes_raw), 0L, df$diabetes_raw))

# Hepamet HFS
df$hfs <- with(df, ifelse(
  !is.na(homa_ir) & !is.na(albumin_gdl),
  calc_hepamet_hfs(sex_male, RIDAGEYR, ast_unl, albumin_gdl, LBXPLTSI, homa_ir, diabetes),
  NA))
df$hfs_high <- as.integer(df$hfs >= 0.47)
cat(sprintf("\n  HFS-eligible (fasting): %d\n", sum(!is.na(df$hfs))))

# ----------------------------------------------------------
# Step 8: FibroScan 子集 (Stack 2, LUX P_ + LUX_J)
# CAP >= 285 (Karlas 2017 + AASLD 2024), LSM >= 8/12 (Castera 2019 + MASLD)
# ----------------------------------------------------------
if ("cap_median" %in% names(df)) {
  df$has_lux <- !is.na(df$cap_median) | !is.na(df$lsm_median)
  df$steatosis_cap <- as.integer(df$cap_median >= 285)
  df$lsm_significant <- as.integer(df$lsm_median >= 8.0)
  df$lsm_advanced <- as.integer(df$lsm_median >= 12.0)
  cat(sprintf("  FibroScan 子集 has_lux: %d / CAP≥285: %d / LSM≥8: %d / LSM≥12: %d\n",
              sum(df$has_lux, na.rm = TRUE),
              sum(df$steatosis_cap, na.rm = TRUE),
              sum(df$lsm_significant, na.rm = TRUE),
              sum(df$lsm_advanced, na.rm = TRUE)))
} else {
  df$has_lux <- FALSE
  cat("  [警告] cap_median 缺 — LUX 子集 N=0\n")
}

# AGILE 3+ + AGILE 4 (FibroScan + age + sex + diabetes + ALT/AST + PLT)
# 公式: AGILE 3+ logit = -19.7 + 0.0454*Age - 0.04*PLT + 1.07*Sex_M + 1.36*Diabetes + 4.61*log10(LSM) - 0.55*AST/ALT
# (简化: Sasso 2022 EASL; 实际 W7 跑时校核公式)
df$agile3_logit <- with(df, ifelse(
  !is.na(lsm_median) & !is.na(LBXPLTSI) & !is.na(ast_unl) & !is.na(alt_unl) & lsm_median > 0,
  -19.7 + 0.0454 * RIDAGEYR - 0.04 * LBXPLTSI + 1.07 * sex_male + 1.36 * diabetes +
    4.61 * log10(lsm_median) - 0.55 * (ast_unl / alt_unl),
  NA))
df$agile3 <- plogis(df$agile3_logit)
df$agile3_advanced <- as.integer(df$agile3 >= 0.451)

# AGILE 4 (cirrhosis prediction, similar formula different coefs)
df$agile4_logit <- with(df, ifelse(
  !is.na(lsm_median) & !is.na(LBXPLTSI) & !is.na(ast_unl) & !is.na(alt_unl) & lsm_median > 0,
  -21.5 + 0.0337 * RIDAGEYR - 0.022 * LBXPLTSI + 0.93 * sex_male + 1.32 * diabetes +
    5.86 * log10(lsm_median) - 0.31 * (ast_unl / alt_unl),
  NA))
df$agile4 <- plogis(df$agile4_logit)
df$agile4_cirrhosis <- as.integer(df$agile4 >= 0.565)

# FAST (FibroScan-AST, Newsome 2020): combine CAP + LSM + AST → NASH + F2+
df$fast_logit <- with(df, ifelse(
  !is.na(cap_median) & !is.na(lsm_median) & !is.na(ast_unl) & lsm_median > 0,
  -1.65 + 1.07 * log(lsm_median) - 2.66e-08 * cap_median^3 + 0.0066 * ast_unl,
  NA))
df$fast <- plogis(df$fast_logit)
df$fast_high <- as.integer(df$fast >= 0.35)
cat(sprintf("  AGILE 3+ advanced (≥0.451): %d / AGILE 4 cirrhosis (≥0.565): %d / FAST high (≥0.35): %d\n",
            sum(df$agile3_advanced, na.rm = TRUE),
            sum(df$agile4_cirrhosis, na.rm = TRUE),
            sum(df$fast_high, na.rm = TRUE)))

# ----------------------------------------------------------
# Step 9: MetALD 三组 (Rinella 2023, 复用 _shared/metald_classify.R)
# Need: steatosis indicator + alcohol g/wk + sex + cardiometabolic risk count
# ----------------------------------------------------------
# 用 FIB-4-advanced 作 steatosis proxy (主); 副用 CAP≥285 for FibroScan 子集
fli_calc <- function(bmi, ggt, tg, waist) {
  ok <- !is.na(bmi) & !is.na(ggt) & !is.na(tg) & !is.na(waist) & tg > 0 & ggt > 0
  L <- ifelse(ok,
              0.953 * log(tg) + 0.139 * bmi + 0.718 * log(ggt) + 0.053 * waist - 15.745,
              NA)
  ifelse(is.na(L), NA, exp(L) / (1 + exp(L)) * 100)
}
df$fli <- fli_calc(df$BMXBMI, df$LBXSGTSI, df$tg_mgdl, df$BMXWAIST)
df$steatosis_fli <- as.integer(!is.na(df$fli) & df$fli >= 60)
# 主 steatosis indicator: 有 LUX 用 CAP, 没 LUX 用 FLI
df$steatosis <- with(df, ifelse(has_lux,
                                steatosis_cap,
                                steatosis_fli))

df$ALQ130 <- na_codes(df$ALQ130, c(77, 99, 777, 999))
df$alcohol_gwk <- alcohol_gwk_from_alq130(df$ALQ130)
df$sex_chr <- ifelse(df$sex_male == 1, "Male", "Female")
df$BPQ020_c <- na_codes(df$BPQ020, c(7, 9))
df$htn_med <- as.integer(!is.na(df$BPQ020_c) & df$BPQ020_c == 1)

df$cm_risk_count <- count_cardiometabolic_risk(
  bmi = df$BMXBMI, waist = df$BMXWAIST, sex = df$sex_chr,
  fbg = df$LBXGLU, hba1c = df$LBXGH,
  t2d_yes = df$diabetes, t2d_med_yes = NA,
  sbp = df$sbp, dbp = df$dbp, htn_med_yes = df$htn_med,
  tg = df$tg_mgdl, lipid_med_yes = NA, hdl = df$hdl_mgdl
)
df$metald_group <- metald_classify(df$steatosis, df$alcohol_gwk, df$sex_chr, df$cm_risk_count)
cat(sprintf("\n  MetALD 三组分布:\n"))
print(table(df$metald_group, useNA = "ifany"))

# ----------------------------------------------------------
# Step 10: Mortality outcome (Stack 1+5)
# ----------------------------------------------------------
df <- df %>% mutate(
  mort_allcause = ifelse(!is.na(MORTSTAT), as.integer(MORTSTAT == 1), NA),
  mort_cm = ifelse(!is.na(MORTSTAT),
                   as.integer(MORTSTAT == 1 & UCOD_LEADING %in% c(1, 5, 7)), NA),
  permth = ifelse(!is.na(PERMTH_EXM), PERMTH_EXM, PERMTH_INT)
)

# ----------------------------------------------------------
# Step 11: 协变量 (COV_M2_COX + 加项)
# ----------------------------------------------------------
df <- df %>% mutate(
  DMDEDUC2 = na_codes(DMDEDUC2, c(7, 9)),
  DMDMARTL = na_codes(DMDMARTL, c(77, 99)),
  SMQ020   = na_codes(SMQ020, c(7, 9)),
  ALQ111   = na_codes(ALQ111, c(7, 9))
)
df <- df %>% mutate(
  age = RIDAGEYR,
  age_group = factor(case_when(
    age >= 20 & age < 40 ~ "20-39",
    age >= 40 & age < 60 ~ "40-59",
    age >= 60            ~ ">=60"
  ), levels = c("20-39","40-59",">=60")),
  race = factor(recode(as.character(RIDRETH1),
                       "1"="Mexican American","2"="Other Hispanic",
                       "3"="Non-Hispanic White","4"="Non-Hispanic Black",
                       "5"="Other Race"),
                levels = c("Non-Hispanic White","Non-Hispanic Black",
                           "Mexican American","Other Hispanic","Other Race")),
  education = factor(case_when(
    DMDEDUC2 %in% c(1,2) ~ "Less than HS",
    DMDEDUC2 == 3        ~ "High school",
    DMDEDUC2 %in% c(4,5) ~ "College or above"
  ), levels = c("Less than HS","High school","College or above")),
  pir = INDFMPIR,
  bmi = BMXBMI,
  smoke = factor(case_when(
    SMQ020 == 2 ~ "Never", SMQ020 == 1 ~ "Ever"
  ), levels = c("Never","Ever")),
  drink = factor(case_when(
    ALQ111 == 1 ~ "Yes", ALQ111 == 2 ~ "No"
  ), levels = c("No","Yes")),
  hypertension = htn_med,
  # 加项: 客观吸烟 cotinine
  cotinine_log = ifelse(!is.na(LBXCOT) & LBXCOT > 0, log(LBXCOT), NA),
  smoke_objective = case_when(
    !is.na(LBXCOT) & LBXCOT >= 10  ~ "Active",       # active smoker
    !is.na(LBXCOT) & LBXCOT >= 0.05 ~ "Passive",     # 二手烟暴露
    !is.na(LBXCOT) & LBXCOT < 0.05 ~ "Non",
    TRUE ~ NA_character_
  )
)
# 加项: Fish freq DRD340-370 (last 30 day fish/shellfish count)
fish_cols <- intersect(c("DRD340","DRD350","DRD360","DRD370"), names(df))
if (length(fish_cols) > 0) {
  for (c in fish_cols) df[[c]] <- na_codes(df[[c]], c(7,9,77,99,777,999))
  df$fish_freq_30d <- rowSums(df[, fish_cols, drop = FALSE], na.rm = TRUE)
}

# 加项: Reproductive 女性 postmenopausal 标志
df$postmenopausal <- with(df, ifelse(
  sex_male == 0 & !is.na(RHQ060) & RHQ060 %in% 1, 1L,
  ifelse(sex_male == 0, 0L, NA_integer_)
))

# 加项: HSCRP inflammation 中介
df$hscrp_log <- ifelse(!is.na(df$LBXHSCRP) & df$LBXHSCRP > 0, log(df$LBXHSCRP), NA)

# ----------------------------------------------------------
# Step 12: Cohort splits + pooled weight (cycle_config.R)
# ----------------------------------------------------------
# 主分析 Stack 1: PFAS D-I + P_ (排除 _J 因 J⊂P_ dedupe; 排除 _C 因 BIOPRO_C 缺)
cycle_years_main <- c(
  NHANES_2005_2006 = 2, NHANES_2007_2008 = 2,
  NHANES_2009_2010 = 2, NHANES_2011_2012 = 2,
  NHANES_2013_2014 = 2, NHANES_2015_2016 = 2,
  PrePandemic_2017_March2020 = 3.2
)

df_main <- df %>% filter(cycle_tag %in% names(cycle_years_main))
log_flow("主分析 Stack 1 cohort (排 _C 无 FIB-4 / 排 _J dedupe)", nrow(df_main))

# Pooled weight — FIX (2026-06-01): PFAS is a ONE-THIRD subsample, so the primary
# analysis weight MUST be the PFAS subsample serum weight (WTSA2YR/WTSC2YR/WTSB2YR/WTSBAPRP),
# NOT the full MEC weight (CDC: "special sample weights are required"). See scripts/_pfas_weight.R.
source("scripts/_pfas_weight.R")
df_main$wt_pooled     <- pooled_pfas_weight(df_main, cycle_years_main)  # PFAS subsample (correct)
df_main$wt_pooled_mec <- pooled_mec_weight(df_main, cycle_years_main)   # retained for reference only
df_main$wt_diet_pooled <- pooled_diet_weight(df_main, cycle_years_main)
df_main$wt_saf_pooled <- pooled_saf_weight(df_main, cycle_years_main)

# Step 13: 核心协变量完整 (drink/cotinine NA 容忍, 用 IPTW/MICE 处理)
core_cov_strict <- c("age","race","education","pir","bmi","RIAGENDR","SDMVPSU","SDMVSTRA")
n_before <- nrow(df_main)
df_main <- df_main %>% filter(if_all(all_of(core_cov_strict), ~ !is.na(.)))
log_flow(sprintf("最终核心协变量完整 (剔 %d)", n_before - nrow(df_main)), nrow(df_main))

# Step 14: Mortality cohort
df_mort <- df_main %>% filter(!is.na(ELIGSTAT), ELIGSTAT == 1)
log_flow("Mortality-linked cohort (Stack 7)", nrow(df_mort))

# Step 15: Stack 2 FibroScan subset (P_ only, has LUX)
df_fibroscan <- df_main %>% filter(has_lux, cycle_tag == "PrePandemic_2017_March2020")
log_flow("Stack 2 FibroScan subset (P_ + has_lux)", nrow(df_fibroscan))

# Step 16: Stack 4 PFAS+Metals subset (PBCD complete)
df_metals <- df_main %>% filter(!is.na(LBXBPB), !is.na(LBXBCD), !is.na(LBXTHG),
                                 !is.na(LBXBSE), !is.na(LBXBMN))
log_flow("Stack 4 PFAS+Metals subset (PBCD 5 metals complete)", nrow(df_metals))

# Metals log + z (同 )
for (m in c("LBXBPB","LBXBCD","LBXTHG","LBXBSE","LBXBMN")) {
  zm <- paste0("z_", tolower(sub("LBX(B|T)", "", m)))
  zm <- sub("HG", "hg", zm); zm <- sub("SE", "se", zm); zm <- sub("MN", "mn", zm)
  df_metals[[paste0("ln_", tolower(sub("LBX(B|T)", "", m)))]] <- log(pmax(df_metals[[m]], 0.01))
  df_metals[[paste0("z_",  tolower(sub("LBX(B|T)", "", m)))]] <- as.numeric(scale(log(pmax(df_metals[[m]], 0.01))))
}

# Save
nhanes_final <- df_main
if (!dir.exists("data/processed")) dir.create("data/processed", recursive = TRUE)
save(nhanes_final, df_mort, df_fibroscan, df_metals,
     cycle_years_main, mort_all, rx_all,
     file = "data/processed/nhanes_final.RData")

cat("\n========================================\n")
cat(sprintf("Stack 1 主分析 N = %d\n", nrow(nhanes_final)))
cat(sprintf("Stack 2 FibroScan N = %d\n", nrow(df_fibroscan)))
cat(sprintf("Stack 4 PFAS+Metals N = %d\n", nrow(df_metals)))
cat(sprintf("Stack 7 Mortality eligible N = %d / 全因死亡 = %d / cardiometabolic 死亡 = %d\n",
            nrow(df_mort),
            sum(df_mort$mort_allcause, na.rm = TRUE),
            sum(df_mort$mort_cm, na.rm = TRUE)))
cat(sprintf("========================================\n"))

if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
flow_df <- do.call(rbind, flow)
write.csv(flow_df, "output/tables/flow_counts.csv", row.names = FALSE)
cat("已保存 data/processed/nhanes_final.RData + output/tables/flow_counts.csv\n")