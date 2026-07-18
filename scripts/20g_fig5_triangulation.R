# ============================================
# 20g_fig5_triangulation.R  (— JHEP Reports 主图 Figure 5)  v2 2026-07-17
#
# 全文主视觉: "同一批人、同一暴露，四个答案"。
#   [1] FIB-4 未调整      -> 六个 PFAS 全部强阳 (指标含 age, 而 PFAS 随龄蓄积)
#   [2] FIB-4 全调整(M2)  -> 全 null
#   [3] APRI 全调整(M2)   -> 四个再度转阳 (指标由 AST 驱动, 无 age 项)
#   [4] 肝硬度 LSM 全调整 -> 全无正相关 (直接测器官, 与两指标无共用项)
# 暴露数据从未变，变的只有测量工具与模型。
#
# v2 相对 v1 的修正:
#   (a) v1 副标题把 M2 误述为 "age, sex and race/ethnicity" —— 那是 M1。M2 是全模型。已修。
#   (b) v1 的 APRI/LSM 95%CI 是从 Wald p 在**百分比尺度**反推 (se=|pct|/z) 的 hack;
#       v2 直接重跑模型取**对数尺度真实 SE**, 再指数变换。CSV 落盘可溯源。
#   (c) 新增 LSM 的 MDES 带 (80% 功效可检出下限), 只画在 LSM 分面 —— 它只适用于 LSM,
#       画满全图会错误暗示 FIB-4/APRI 也是这个功效。
#
# 输出: output/tables/table_triangulation.csv  (图的唯一数据源, 全精度)
#       output/figures/fig5_triangulation.png  (300 dpi)
# ============================================
suppressPackageStartupMessages({ library(survey); library(dplyr); library(ggplot2) })
options(survey.lonely.psu = "adjust")
load("data/processed/nhanes_final.RData"); load("data/processed/nhanes_design.RData")

pf6 <- c("LBXPFNA_z","LBXPFOA_z","LBXPFHS_z","LBXPFOS_z","LBXPFDE_z","LBXMPAH_z")
nm  <- c(LBXPFNA_z="PFNA", LBXPFOA_z="PFOA", LBXPFHS_z="PFHxS",
         LBXPFOS_z="PFOS", LBXPFDE_z="PFDA", LBXMPAH_z="Me-PFOSA-AcOH")
ord <- unname(nm[pf6])
M2  <- c("age","RIAGENDR","race","education","pir","smoke","bmi","diabetes","hypertension")

# 取对数尺度 estimate + SE，再统一指数变换成 %/SD
grab <- function(des, outcome, pf, covs) {
  f <- as.formula(paste(outcome, "~", paste(c(pf, covs), collapse = " + ")))
  m <- tryCatch(svyglm(f, design = des), error = function(e) NULL)
  if (is.null(m) || !(pf %in% rownames(summary(m)$coef))) return(NULL)
  s  <- summary(m)$coef[pf, ]
  b  <- unname(s["Estimate"]); se <- unname(s["Std. Error"])
  data.frame(pfas = nm[[pf]], b = b, se = se, p = unname(s["Pr(>|t|)"]),
             n = nobs(m), stringsAsFactors = FALSE)
}

# --- 三个测量各自的设计 ---
des_main <- update(design_main, log_apri = log(pmax(apri, 0.01)))
fv  <- design_fibroscan$variables
des_lux <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTSBAPRP, data = fv, nest = TRUE)
des_lux <- update(des_lux, log_lsm = log(pmax(lsm_median, 0.1)))

cvM2_main <- intersect(M2, names(des_main$variables))
cvM2_lux  <- intersect(M2, names(des_lux$variables))
stopifnot("age 必须在协变量中" = "age" %in% cvM2_main, "RIAGENDR 必须在" = "RIAGENDR" %in% cvM2_main)
cat("M2 协变量 (main):", paste(cvM2_main, collapse=", "), "\n")
cat("M2 协变量 (lux) :", paste(cvM2_lux,  collapse=", "), "\n\n")

rows <- list()
for (pf in pf6) {
  rows[[length(rows)+1]] <- cbind(grab(des_main, "fib4_log", pf, character(0)),
                                  measure = "FIB-4", model = "Unadjusted")
  rows[[length(rows)+1]] <- cbind(grab(des_main, "fib4_log", pf, cvM2_main),
                                  measure = "FIB-4", model = "Fully adjusted")
  rows[[length(rows)+1]] <- cbind(grab(des_main, "log_apri", pf, cvM2_main),
                                  measure = "APRI", model = "Fully adjusted")
  rows[[length(rows)+1]] <- cbind(grab(des_lux,  "log_lsm",  pf, cvM2_lux),
                                  measure = "Liver stiffness", model = "Fully adjusted")
}
dat <- bind_rows(rows) %>%
  mutate(pct = (exp(b) - 1) * 100,
         lo  = (exp(b - 1.96 * se) - 1) * 100,
         hi  = (exp(b + 1.96 * se) - 1) * 100,
         mdes_pct = (exp((qnorm(0.975) + qnorm(0.80)) * se) - 1) * 100)

write.csv(dat, "output/tables/table_triangulation.csv", row.names = FALSE)
cat("[OK] output/tables/table_triangulation.csv\n")

# LSM 的 MDES 范围 —— 只用于 LSM 分面的功效带
mdes_rng <- dat %>% filter(measure == "Liver stiffness") %>% summarise(lo = min(mdes_pct), hi = max(mdes_pct))
cat(sprintf("LSM MDES 范围: %.1f%% – %.1f%% /SD\n", mdes_rng$lo, mdes_rng$hi))

panel_lab <- c(
  "FIB-4\nUnadjusted"      = "1. FIB-4, unadjusted",
  "FIB-4\nFully adjusted"  = "2. FIB-4, fully adjusted",
  "APRI\nFully adjusted"   = "3. APRI, fully adjusted",
  "Liver stiffness\nFully adjusted" = "4. Liver stiffness, fully adjusted")
dat <- dat %>%
  mutate(panel = factor(paste0(measure, "\n", model),
                        levels = names(panel_lab), labels = unname(panel_lab)),
         pfas  = factor(pfas, levels = rev(ord)),
         sig   = ifelse(p < 0.05, "P < 0.05", "Not significant"))

# 不画 MDES 功效带: ±4.9% 超出第4面板可视范围, 视觉上像"整面板不可分辨"(夸大弱点);
# 且 MDES 是 80% 功效下限而非硬阈值, 带内仍可能显著(PFHxS p=0.017 就在带内) → 图上说不清。
# 功效已在摘要/Results §3.6/Discussion/Table S9 四处交代, 图注给数字即可。
p5 <- ggplot(dat, aes(x = pct, y = pfas)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi, colour = panel), height = 0.20, linewidth = 0.55) +
  geom_point(aes(colour = panel, shape = sig), size = 2.4, fill = "white", stroke = 0.8) +
  facet_grid(. ~ panel, scales = "free_x") +
  scale_shape_manual(values = c("P < 0.05" = 16, "Not significant" = 21), name = NULL) +
  scale_colour_manual(values = setNames(c("#B03A2E", "#7F8C8D", "#E67E22", "#2471A3"), unname(panel_lab)),
                      guide = "none") +
  labs(
    x = "Change in fibrosis measure per SD increment in PFAS (%)",
    y = NULL,
    title = "Same adults, same exposures, four answers",
    subtitle = paste0("Panels 2–4 apply one covariate set (age, sex, race/ethnicity, education, income, smoking, BMI, ",
                      "diabetes, hypertension) and differ\nin how fibrosis was measured. Panel 4 is restricted to the ",
                      "1,929 of these adults who underwent elastography."),
    caption = paste0(
      "n = 8,588 (7,454 for PFOS and PFOA, which exclude the 2013–2014 cycle) in panels 2–3; 8,590 and 7,456 unadjusted; ",
      "1,929 in panel 4.\n",
      "Panel 1 has its own x-axis scale; panels 2–4 share one. Filled circles P < 0.05; open circles not significant. ",
      "95% CIs from model standard errors on the log scale.\n",
      "Panel 4 resolves associations of ", sprintf("%.1f–%.1f%%", mdes_rng$lo, mdes_rng$hi),
      " per SD at 80% power (Table S9), so it withholds corroboration from the panel-3 signal rather than excluding it.")
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.spacing = unit(0.9, "lines"),
        strip.text = element_text(face = "bold", size = 9.6),
        legend.position = "top",
        plot.title = element_text(face = "bold", size = 13.5),
        plot.subtitle = element_text(size = 9, colour = "grey35"),
        plot.caption = element_text(size = 7.2, colour = "grey40", hjust = 0),
        axis.text.y = element_text(size = 9.6))

ggsave("output/figures/fig5_triangulation.png", plot = p5,
       width = 11.0, height = 4.6, dpi = 300, bg = "white")
cat("[OK] output/figures/fig5_triangulation.png\n\n")
print(dat %>% select(pfas, panel, pct, lo, hi, p, n) %>% arrange(panel, pfas), digits = 3)