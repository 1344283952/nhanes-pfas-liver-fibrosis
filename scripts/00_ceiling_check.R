# ============================================
# scripts/00_ceiling_check.R
# Phase 1a: NHANES PFAS × MASLD/FibroScan + Mortality 子集天花板实查
#
# 任务: 在写 task.md 前用真实 R 实查 N (按 [feedback-nhanes-sample-ceiling-first] 教训)
# 输出: output/_ceiling_report.csv (7 维 stack N)
# 资源: 10 worker 并行下载 ~50 个 .xpt (~3-5 min); 不抢 /BKMR (单 R 进程)
# 触发: cd  && Rscript scripts/00_ceiling_check.R
# ============================================

cat("==========================================\n")
cat("Phase 1a: PFAS 子集天花板实查 ()\n")
cat("==========================================\n\n")

suppressPackageStartupMessages({
  library(parallel)
  library(haven)
  library(dplyr)
})

# 输出 + 临时目录
tmp_dir <- "data/raw/_ceiling_tmp"
if (!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = TRUE)
if (!dir.exists("output")) dir.create("output", recursive = TRUE)

# ----------------------------------------------------------
# 下载 helper (内容嗅探, 复用 templates/01_download_data.R 风格)
# ----------------------------------------------------------
download_nhanes <- function(year, fname, dest_dir) {
  fp <- file.path(dest_dir, paste0(fname, ".xpt"))
  is_real_xpt <- function(p) {
    if (!file.exists(p) || file.size(p) < 1024) return(FALSE)
    con <- file(p, "rb"); on.exit(close(con))
    identical(readChar(con, 6, useBytes = TRUE), "HEADER")
  }
  if (is_real_xpt(fp)) return(list(file = fname, ok = TRUE))
  for (u in c(
    sprintf("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/%d/DataFiles/%s.xpt", year, fname),
    sprintf("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/%d/DataFiles/%s.XPT", year, fname)
  )) {
    ok <- tryCatch({
      download.file(u, fp, mode = "wb", quiet = TRUE, method = "libcurl")
      is_real_xpt(fp)
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (isTRUE(ok)) return(list(file = fname, ok = TRUE))
    if (file.exists(fp)) file.remove(fp)
  }
  list(file = fname, ok = FALSE)
}

# ----------------------------------------------------------
# 任务清单
# PFAS: C(2003) D(2005) E(2007) F(2009) G(2011) H(2013) I(2015) J(2017) + P_(2017-Mar2020)
# Cross: DEMO/BIOPRO/CBC/HSCRP/PBCD/IHGEM/UHM/PHTHTE/DXX/DR1TOT/RHQ/COT/HEPB_S/HEPC/LUX/ALB_CR
# ----------------------------------------------------------
pfas_modules <- data.frame(
  label = c("PFAS_C","PFAS_D","PFAS_E","PFAS_F","PFAS_G","PFAS_H","PFAS_I","PFAS_J","PFAS_PP"),
  year  = c(2003, 2005, 2007, 2009, 2011, 2013, 2015, 2017, 2017),
  file  = c("L24PFC_C","PFC_D","PFC_E","PFC_F","PFC_G","PFAS_H","PFAS_I","PFAS_J","P_PFAS"),
  stringsAsFactors = FALSE
)

cycle_suffix <- c(C=2003, D=2005, E=2007, F=2009, G=2011, H=2013, I=2015, J=2017)
modules_cross <- c("DEMO","BIOPRO","CBC","HSCRP","PBCD","IHGEM","UHM",
                   "PHTHTE","DXX","DR1TOT","RHQ","COT","HEPB_S","HEPC","ALB_CR")
cross_tasks <- do.call(rbind, lapply(modules_cross, function(m) {
  data.frame(label = paste0(m, "_", names(cycle_suffix)),
             year = unname(cycle_suffix),
             file = paste0(m, "_", names(cycle_suffix)),
             stringsAsFactors = FALSE)
}))

# Pre-pandemic P_ 文件
pp_modules <- c("DEMO","BIOPRO","CBC","HSCRP","PBCD","PHTHTE","DXX",
                "DR1TOT","RHQ","COT","LUX","ALB_CR")
pp_tasks <- data.frame(
  label = paste0(pp_modules, "_PP"),
  year = 2017,
  file = paste0("P_", pp_modules),
  stringsAsFactors = FALSE
)

# LUX (FibroScan) 仅 J + P_ (非全周期)
lux_tasks <- data.frame(label = "LUX_J", year = 2017, file = "LUX_J", stringsAsFactors = FALSE)

all_tasks <- rbind(pfas_modules, cross_tasks, pp_tasks, lux_tasks)
n_tasks <- nrow(all_tasks)
cat(sprintf("总下载任务: %d (PFAS=%d, cross-cycle=%d, P_=%d, LUX_J=%d)\n",
            n_tasks, nrow(pfas_modules), nrow(cross_tasks), nrow(pp_tasks), nrow(lux_tasks)))

# ----------------------------------------------------------
# 并行下载 (PSOCK 10 worker, ~3-5 min)
# ----------------------------------------------------------
cat("\n启动 10 worker 并行下载...\n")
# 把每个 task 完整 list 传给 worker, 不依赖主进程 all_tasks (PSOCK 隔离)
task_list <- lapply(seq_len(n_tasks), function(i) {
  list(year = all_tasks$year[i], file = all_tasks$file[i])
})
cl <- makeCluster(10)
clusterExport(cl, c("download_nhanes", "tmp_dir"), envir = environment())
t0 <- Sys.time()
results <- parLapplyLB(cl, task_list, function(t) {
  download_nhanes(t$year, t$file, tmp_dir)
})
stopCluster(cl)
elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
all_tasks$downloaded <- vapply(results, function(r) isTRUE(r$ok), logical(1))
cat(sprintf("下载完成 %d/%d (%.2f min)\n", sum(all_tasks$downloaded), n_tasks, elapsed))

fails <- all_tasks$file[!all_tasks$downloaded]
if (length(fails) > 0) {
  cat("失败 (部分模块周期不存在是正常的, 例如 PFAS_C 历史命名差异):\n")
  for (f in fails) cat("  -", f, "\n")
}

# ----------------------------------------------------------
# Helpers: 读取 + 求 PFAS 有效 SEQN
# ----------------------------------------------------------
load_xpt_safe <- function(fname) {
  fp <- file.path(tmp_dir, paste0(fname, ".xpt"))
  if (!file.exists(fp)) return(NULL)
  tryCatch(haven::read_xpt(fp), error = function(e) NULL)
}

pfas_seqn <- function(pfas_df) {
  if (is.null(pfas_df)) return(integer(0))
  pf_cols <- grep("^(LBX|LBD)PF", names(pfas_df), value = TRUE)
  if (length(pf_cols) == 0) return(integer(0))
  pfas_df$SEQN[rowSums(!is.na(pfas_df[, pf_cols, drop = FALSE])) > 0]
}

# 跨周期 PFAS file → DEMO/BIOPRO/CBC file 映射
file_map <- list(
  PFAS_C = list(demo="DEMO_C", biopro="BIOPRO_C", cbc="CBC_C"),
  PFAS_D = list(demo="DEMO_D", biopro="BIOPRO_D", cbc="CBC_D"),
  PFAS_E = list(demo="DEMO_E", biopro="BIOPRO_E", cbc="CBC_E"),
  PFAS_F = list(demo="DEMO_F", biopro="BIOPRO_F", cbc="CBC_F"),
  PFAS_G = list(demo="DEMO_G", biopro="BIOPRO_G", cbc="CBC_G"),
  PFAS_H = list(demo="DEMO_H", biopro="BIOPRO_H", cbc="CBC_H"),
  PFAS_I = list(demo="DEMO_I", biopro="BIOPRO_I", cbc="CBC_I"),
  PFAS_J = list(demo="DEMO_J", biopro="BIOPRO_J", cbc="CBC_J"),
  PFAS_PP = list(demo="P_DEMO", biopro="P_BIOPRO", cbc="P_CBC")
)

# ----------------------------------------------------------
# Stack 0+1: PFAS × DEMO ≥20 / PFAS × FIB-4(ALT/AST/PLT) ≥20  按周期
# ----------------------------------------------------------
cat("\n求 7 维 ceiling 交集 N...\n\n")

results_by_cycle <- lapply(pfas_modules$label, function(lab) {
  pf <- load_xpt_safe(pfas_modules$file[pfas_modules$label == lab])
  fm <- file_map[[lab]]
  demo <- load_xpt_safe(fm$demo); biopro <- load_xpt_safe(fm$biopro); cbc <- load_xpt_safe(fm$cbc)
  ps <- pfas_seqn(pf)
  n_pfas <- length(ps)
  n_pfas_adult <- if (!is.null(demo)) sum(demo$SEQN %in% ps & demo$RIDAGEYR >= 20) else NA
  n_fib4 <- NA
  if (!is.null(demo) && !is.null(biopro) && !is.null(cbc)) {
    alt_col <- intersect(c("LBXSATSI","LBXSALT"), names(biopro))[1]
    ast_col <- intersect(c("LBXSASSI","LBXSAST"), names(biopro))[1]
    if (!is.na(alt_col) && !is.na(ast_col) && "LBXPLTSI" %in% names(cbc)) {
      lf_seqn <- biopro$SEQN[!is.na(biopro[[alt_col]]) & !is.na(biopro[[ast_col]])]
      pl_seqn <- cbc$SEQN[!is.na(cbc$LBXPLTSI)]
      inter <- Reduce(intersect, list(ps, lf_seqn, pl_seqn))
      n_fib4 <- sum(demo$SEQN %in% inter & demo$RIDAGEYR >= 20)
    }
  }
  data.frame(cycle = lab, n_pfas_total = n_pfas,
             n_pfas_adult = n_pfas_adult, n_fib4_adult = n_fib4,
             stringsAsFactors = FALSE)
})
stack_01 <- do.call(rbind, results_by_cycle)
cat("Stack 0 + 1 — PFAS / PFAS×FIB-4 各周期 ≥20 岁 N:\n")
print(stack_01)

n_fib4_pooled <- sum(stack_01$n_fib4_adult[stack_01$cycle != "PFAS_J"], na.rm = TRUE)
cat(sprintf("\n→ PFAS C-J + P_ pooled FIB-4 (J⊂P_ dedupe 排除 PFAS_J 行): %d\n", n_fib4_pooled))

# ----------------------------------------------------------
# Stack 2 — PFAS × FibroScan (LUX P_ only)
# ----------------------------------------------------------
pf_pp <- load_xpt_safe("P_PFAS")
lux_pp <- load_xpt_safe("P_LUX")
demo_pp <- load_xpt_safe("P_DEMO")
n_stack2 <- NA
if (!is.null(pf_pp) && !is.null(lux_pp) && !is.null(demo_pp)) {
  ps_pp <- pfas_seqn(pf_pp)
  lux_valid <- lux_pp$SEQN[!is.na(lux_pp$LUXSMED) | !is.na(lux_pp$LUXCAPM)]
  n_stack2 <- sum(demo_pp$SEQN %in% Reduce(intersect, list(ps_pp, lux_valid)) & demo_pp$RIDAGEYR >= 20)
}
cat(sprintf("Stack 2 — PFAS_P × LUX (LSM/CAP non-NA) × ≥20: %s\n",
            if (!is.na(n_stack2)) n_stack2 else "P_LUX 缺失"))

# ----------------------------------------------------------
# Stack 3 — PFAS × Metals (PBCD: Pb)
# ----------------------------------------------------------
pbcd_pp <- load_xpt_safe("P_PBCD")
n_stack3 <- NA
if (!is.null(pf_pp) && !is.null(pbcd_pp) && !is.null(demo_pp)) {
  metal_seqn <- pbcd_pp$SEQN[!is.na(pbcd_pp$LBXBPB)]
  ps_pp <- pfas_seqn(pf_pp)
  n_stack3 <- sum(demo_pp$SEQN %in% Reduce(intersect, list(ps_pp, metal_seqn)) & demo_pp$RIDAGEYR >= 20)
}
cat(sprintf("Stack 3 — PFAS_P × PBCD (Pb non-NA) × ≥20: %s\n",
            if (!is.na(n_stack3)) n_stack3 else "P_PBCD 缺失"))

# ----------------------------------------------------------
# Stack 4 — PFAS × Phthalates (P_PHTHTE 存在性 关键 Agent F 不确定)
# ----------------------------------------------------------
phth_pp <- load_xpt_safe("P_PHTHTE")
n_stack4 <- NA
if (!is.null(pf_pp) && !is.null(phth_pp) && !is.null(demo_pp)) {
  phth_cols <- grep("^URX", names(phth_pp), value = TRUE)
  if (length(phth_cols) > 0) {
    phth_seqn <- phth_pp$SEQN[rowSums(!is.na(phth_pp[, phth_cols, drop = FALSE])) > 0]
    ps_pp <- pfas_seqn(pf_pp)
    n_stack4 <- sum(demo_pp$SEQN %in% Reduce(intersect, list(ps_pp, phth_seqn)) & demo_pp$RIDAGEYR >= 20)
  }
}
cat(sprintf("Stack 4 — PFAS_P × Phthalates × ≥20: %s  [Agent F 不确定项, 实查]\n",
            if (!is.na(n_stack4)) n_stack4 else "P_PHTHTE 缺失或无数据"))

# ----------------------------------------------------------
# Stack 5 — PFAS × DXA (P_DXX 存在性 关键 Agent F 不确定)
# ----------------------------------------------------------
dxx_pp <- load_xpt_safe("P_DXX")
n_stack5 <- NA
if (!is.null(pf_pp) && !is.null(dxx_pp) && !is.null(demo_pp)) {
  dxa_cols <- grep("^DX", names(dxx_pp), value = TRUE)
  if (length(dxa_cols) > 0) {
    dxa_seqn <- dxx_pp$SEQN[rowSums(!is.na(dxx_pp[, dxa_cols, drop = FALSE])) > 0]
    ps_pp <- pfas_seqn(pf_pp)
    n_stack5 <- sum(demo_pp$SEQN %in% Reduce(intersect, list(ps_pp, dxa_seqn)) & demo_pp$RIDAGEYR >= 20)
  }
}
cat(sprintf("Stack 5 — PFAS_P × DXA × ≥20: %s  [Agent F 不确定项, 实查]\n",
            if (!is.na(n_stack5)) n_stack5 else "P_DXX 缺失"))

# ----------------------------------------------------------
# Stack 6 — PFAS_J × HEPB/HEPC 排除 (P_ cycle HEPB/HEPC 文件状态待 W2)
# ----------------------------------------------------------
hep_b <- load_xpt_safe("HEPB_S_J"); hep_c <- load_xpt_safe("HEPC_J")
pfas_j <- load_xpt_safe("PFAS_J")
n_stack6 <- NA
if (!is.null(hep_b) && !is.null(hep_c) && !is.null(pfas_j)) {
  hbs_pos <- hep_b$SEQN[!is.na(hep_b$LBXHBS) & hep_b$LBXHBS == 1]
  hcv_pos <- hep_c$SEQN[!is.na(hep_c$LBXHCR) & hep_c$LBXHCR == 1]
  ps_j <- pfas_seqn(pfas_j)
  n_stack6 <- length(intersect(ps_j, union(hbs_pos, hcv_pos)))
}
cat(sprintf("Stack 6 — PFAS_J × HEPB/HEPC 阳性 (需排除): %s\n",
            if (!is.na(n_stack6)) n_stack6 else "HEPB_S_J/HEPC_J 缺失"))

# ----------------------------------------------------------
# Stack 7 — PFAS × Mortality LMF (NDI 2019, 假设 99% link 按 实测)
# ----------------------------------------------------------
n_pfas_mort_eligible <- round(n_fib4_pooled * 0.99)
cat(sprintf("Stack 7 — PFAS C-J × LMF 2019 eligible (~99%% 按 ): %d\n", n_pfas_mort_eligible))
cat("  估算事件 (按 全因 ~6%% / cardiometabolic ~2%%):\n")
cat(sprintf("    全因死亡: ~%d\n", round(n_pfas_mort_eligible * 0.06)))
cat(sprintf("    Cardiometabolic 死亡: ~%d\n", round(n_pfas_mort_eligible * 0.02)))

# ----------------------------------------------------------
# 输出 _ceiling_report.csv
# ----------------------------------------------------------
ceiling_report <- rbind(
  stack_01 |> mutate(stack_type = "PFAS_FIB4_BY_CYCLE") |>
    select(stack_type, stack = cycle, n_pfas_total, n_pfas_adult, n_intersect = n_fib4_adult),
  data.frame(
    stack_type = c("PFAS_FIB4_POOLED_DEDUPE", "STACK2_FIBROSCAN_P", "STACK3_METALS_P",
                   "STACK4_PHTHALATES_P", "STACK5_DXA_P", "STACK6_HEP_EXCLUSION_J",
                   "STACK7_MORTALITY_ELIGIBLE"),
    stack = c("C-J+P_ dedupe", "PFAS_P×LUX×≥20", "PFAS_P×PBCD×≥20",
              "PFAS_P×PHTHTE×≥20", "PFAS_P×DXX×≥20", "PFAS_J×HEPpos",
              "PFAS_C-J pooled × LMF 99%"),
    n_pfas_total = NA, n_pfas_adult = NA,
    n_intersect = c(n_fib4_pooled, n_stack2, n_stack3, n_stack4, n_stack5,
                    n_stack6, n_pfas_mort_eligible),
    stringsAsFactors = FALSE
  )
)
write.csv(ceiling_report, "output/_ceiling_report.csv", row.names = FALSE)
cat("\n==========================================\n")
cat("✓ _ceiling_report.csv 已写入 output/\n")
cat("✓ 把上面 7 个 Stack 数字交给 task.md v1 (Phase 1b)\n")
cat("==========================================\n")