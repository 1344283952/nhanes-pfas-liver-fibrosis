# ============================================
# 06e_ph_check.R — 比例风险(PH)假定检验
#
# 为什么单独一个脚本 + 为什么 unweighted:
#   cox.zph() 作用在 svycoxph 对象上会产生**退化结果**（chisq≈1e-6, p≈1.0），
#   看起来"完全不违反 PH"，实为假象——加权 Cox 的 Schoenfeld 残差尺度被权重扭曲。
#   正确做法: 用 **unweighted complete-case coxph** 做 PH 诊断。
#   PH 检验是对"模型设定"的诊断，不是对总体的推断，因此不需要抽样权重；
#   加权估计仍以 06d 的 svycoxph 为准，本脚本只回答"PH 是否成立"。
#   (已废弃的 06_table2.R:159 正是 cox.zph(svycoxph(...)) 的错误写法，勿复用。)
#
# 输入: data/processed/nhanes_final.RData, nhanes_design.RData（与 06d 同源）
# 输出: output/tables/cox_schoenfeld_pfas.csv
# ============================================
suppressPackageStartupMessages({ library(survey); library(survival); library(dplyr) })
options(survey.lonely.psu = "adjust")
load("data/processed/nhanes_final.RData"); load("data/processed/nhanes_design.RData")

clean_nm <- c(LBXPFNA="PFNA", LBXMPAH="Me-PFOSA-AcOH", LBXPFDE="PFDA",
              LBXPFOS="PFOS", LBXPFOA="PFOA", LBXPFHS="PFHxS")
pf6 <- c("LBXPFNA_z","LBXMPAH_z","LBXPFDE_z","LBXPFOS_z","LBXPFOA_z","LBXPFHS_z")

# 与 06d 完全一致的 M2（confounders only，匹配 DAG 与正文）
M2_clean <- c("age","RIAGENDR","race","education","pir","smoke","bmi","diabetes","hypertension")
keepcov <- function(df, covs) covs[vapply(covs, function(c) {
  x <- df[[c]]; mean(is.na(x)) <= 0.20 && length(unique(x[!is.na(x)])) >= 2 }, logical(1))]

desm <- design_mortality
desm$variables$time_yr <- desm$variables$permth / 12
vm <- desm$variables
covM2m <- keepcov(vm, intersect(M2_clean, names(vm)))

# --- 护栏: 核对协变量真进了模型（教训: keepcov 曾把空 sex 悄悄丢掉 → 从未校正性别）---
cat("PH 检验使用的 M2 协变量:", paste(covM2m, collapse=", "), "\n")
dropped <- setdiff(M2_clean, covM2m)
if (length(dropped)) cat("!! 被 keepcov 丢弃:", paste(dropped, collapse=", "), "\n")
stopifnot("age 必须在协变量中" = "age" %in% covM2m)
stopifnot("RIAGENDR(性别) 必须在协变量中" = "RIAGENDR" %in% covM2m)

cat(sprintf("死亡队列 n=%d | 全因死亡=%d | 心代谢死亡=%d\n\n",
            nrow(vm), sum(vm$mort_allcause == 1, na.rm = TRUE), sum(vm$mort_cm == 1, na.rm = TRUE)))

rows <- list()
for (pf in pf6) {
  for (oc in c("mort_allcause", "mort_cm")) {
    f <- as.formula(sprintf("Surv(time_yr, %s) ~ %s + %s", oc, pf, paste(covM2m, collapse = " + ")))
    # unweighted complete-case —— 关键: 不用 svycoxph
    dat <- vm[, c("time_yr", oc, pf, covM2m)]
    dat <- dat[complete.cases(dat), ]
    fit <- tryCatch(coxph(f, data = dat), error = function(e) NULL)
    if (is.null(fit)) { cat("FAIL:", pf, oc, "\n"); next }
    zph <- tryCatch(cox.zph(fit), error = function(e) NULL)
    if (is.null(zph)) { cat("zph FAIL:", pf, oc, "\n"); next }
    tb <- zph$table
    rows[[length(rows) + 1]] <- data.frame(
      pfas          = clean_nm[[sub("_z$", "", pf)]],
      outcome       = oc,
      n_completecase= nrow(dat),
      events        = sum(dat[[oc]] == 1, na.rm = TRUE),
      term          = "PFAS",
      chisq_pfas    = round(tb[pf, "chisq"], 4),
      p_pfas        = round(tb[pf, "p"], 4),
      chisq_global  = round(tb["GLOBAL", "chisq"], 4),
      p_global      = round(tb["GLOBAL", "p"], 4),
      stringsAsFactors = FALSE)
  }
}
out <- do.call(rbind, rows)
write.csv(out, "output/tables/cox_schoenfeld_pfas.csv", row.names = FALSE)
cat("=== Schoenfeld PH 检验 (unweighted complete-case) ===\n"); print(out)

# --- 退化检测: 若 chisq 全部 ~0 且 p 全部 ~1，说明又跑成了加权版的假象 ---
if (all(out$chisq_pfas < 1e-4) && all(out$p_pfas > 0.99)) {
  stop("PH 结果退化 (chisq≈0, p≈1) —— 检查是否误用了加权模型")
}
viol <- out[out$p_pfas < 0.05, ]
cat(sprintf("\nPFAS 项违反 PH 的模型数: %d / %d\n", nrow(viol), nrow(out)))
if (nrow(viol)) { cat("违反 PH 的:\n"); print(viol[, c("pfas","outcome","p_pfas")]) }
cat(sprintf("整体(GLOBAL)违反 PH 的模型数: %d / %d\n", sum(out$p_global < 0.05), nrow(out)))
cat("\n输出 -> output/tables/cox_schoenfeld_pfas.csv\n")