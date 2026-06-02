# ============================================
# run_all.R — reproduce the reported analysis (compound-resolved PFAS x liver fibrosis)
# 用法: 在 projects/007.../ 目录下执行  Rscript scripts/run_all.R
# 若已有 data/processed/nhanes_final.RData + nhanes_design.RData, 可从步骤 5 (06d) 开始。
# ============================================

cat("== NHANES PFAS x liver fibrosis — full reproduction ==\n\n")
start_time <- Sys.time()

# 完整性预检: AST/ALT 映射 + FIB-4 公式方向 —— 映射错就 stop(), 不会跑出假结果
# (2026-06-02 新增, 防 那种"FIB-4 整反"的静默致命 bug 复发; 见 _integrity_preflight.R)
.ipf <- "templates/_shared/_integrity_preflight.R"
if (file.exists(.ipf)) source(.ipf) else cat("[warn] 未找到 _integrity_preflight.R, 跳过完整性预检\n")

scripts <- c(
  "scripts/01_download_data.R",         # 1. download NHANES PFAS + biochemistry + mortality
  "scripts/02_merge_data.R",            # 2. merge multi-cycle (2013-14 PFOS/PFOA excluded by design)
  "scripts/03_clean_data.R",            # 3. clean + isomer reconstruction + PFAS-subsample weight
  "scripts/04_survey_design.R",         # 4. complex-survey design (PFAS-subsample serum weights)
  "scripts/06d_reanalysis_cleanM2.R",   # 5. AUTHORITATIVE: single/mutual/threshold + mortality + E-value + BH(4 families) + distribution
  "scripts/07_table3.R",                # 6. mixture (qgcomp + WQS)
  "scripts/_agefree_check.R",           # 6b. age-free APRI triangulation (Table S7)
  "scripts/_fibroscan_check.R",         # 6c. FibroScan liver-stiffness triangulation (Table S7)
  "scripts/20e_figures_corrected.R",    # 7. Figures 1-4 + Fig S1 (read from the CSVs above)
  "scripts/20f_graphabs_table1.R",      # 8. graphical abstract (+ short Table 1, superseded next)
  "scripts/_table1_full.R",             # 9. expanded Table 1 by PFNA tertile (final table1_by_pfna.csv)
  "scripts/_render_citations.R"         # 10. render Vancouver citations into manuscript_v2.md
)

for (s in scripts) {
  cat(paste0("\n>>> ", s, "\n"))
  tryCatch({ source(s, echo = FALSE); cat(paste0(">>> ", s, " OK\n")) },
           error = function(e) cat(paste0(">>> ", s, " FAILED: ", e$message, "\n")))
}

elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 1)
cat(sprintf("\n== Done in %s min ==\n", elapsed))
cat("Key outputs: output/tables/table_xsec_pfas_fib4.csv (single-PFAS), table_mutual_pfas_fib4.csv (joint),\n")
cat("  evalue.csv, table_bh_fdr.csv, table2_cox_pfas.csv, fine_gray_cm_pfas.csv, table1_by_pfna.csv,\n")
cat("  mixture_triangulation_stack1.csv; figures fig1_consort..fig4_threshold_or + figS_pfas_corr + graphical_abstract.\n")