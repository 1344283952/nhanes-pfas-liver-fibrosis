# ============================================
# scripts/01_download_data.R
# 下载 NHANES PFAS C-J + P_PFAS (2003-March 2020) + 12 cross 模块 + NCHS LMF
#
# 论文: PFAS × MASLD/FibroScan + Cardiometabolic Mortality
#       主投 J Hepatol IF 26.8
# Cohort:
#   Stack 1 主: PFAS C-I + P_ (J⊂P_ dedupe), age ≥ 20 + FIB-4, N=12,533
#   Stack 2 副: PFAS_P × LUX (P_ only), N=2,391
#   Stack 4 联合: PFAS+Metals (PBCD/IHGEM/UHM), 复用 框架
#
# 8 worker 并行 (CDC 限 ~10), ~3-5 min 完成
# ============================================

cat("========================================\n")
cat("NHANES PFAS C-J + P_PFAS + 12 模块下载\n")
cat("========================================\n\n")

raw_dir  <- "data/raw"
mort_dir <- "data/raw/mortality"
log_path <- "data/raw/_download.log"
if (!dir.exists(raw_dir))  dir.create(raw_dir,  recursive = TRUE)
if (!dir.exists(mort_dir)) dir.create(mort_dir, recursive = TRUE)

writeLines(c("# NHANES download log", paste("# started:", Sys.time())),
           log_path)
log_line <- function(msg) {
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = log_path, append = TRUE)
}

# --------------------------------------------------
# 周期: C(2003) D(2005) E(2007) F(2009) G(2011) H(2013) I(2015) J(2017) + P_
# --------------------------------------------------
cycles_std <- data.frame(
  suffix   = c("_C","_D","_E","_F","_G","_H","_I","_J"),
  pathYear = c( 2003, 2005, 2007, 2009, 2011, 2013, 2015, 2017),
  stringsAsFactors = FALSE
)

modules_uniform <- c(
  # 核心人口学 + 临床
  "DEMO","BMX","BPX","SMQ","ALQ","DIQ","BPQ","MCQ","PAQ","RHQ",
  "RXQ_RX",  # 长表
  # 暴露 (Stack 4 mixture + Stack 1 协变)
  "PBCD","IHGEM","UHM",
  # 结局 (Stack 2 FibroScan; LUX 仅 J cycle 有, 但加 _ 各 cycle 试 → 缺的自动失败)
  "LUX",
  # FIB-4 + 5-score (AGILE + FAST + Hepamet)
  "BIOPRO","CBC","HSCRP",
  # HOMA-IR
  "GHB","GLU","INS",
  # 血脂 (5-score + cardiometabolic)
  "HDL","TCHOL","TRIGLY",
  # 膳食 (fish freq DRD340-370 + DII covariate)
  "DR1TOT","DR2TOT",
  # 协变量
  "VID","ALB_CR","COT",
  # viral hepatitis 排除
  "HEPB_S","HEPC"
)

# PFAS 跨周期文件名漂移 (CDC 命名不一致)
pfas_files <- data.frame(
  filename = c("L24PFC_C","PFC_D","PFC_E","PFC_F","PFC_G",
               "PFAS_H","PFAS_I","PFAS_J","P_PFAS"),
  pathYear = c( 2003,     2005,   2007,   2009,   2011,
                2013,     2015,   2017,    2017),
  stringsAsFactors = FALSE
)

# --------------------------------------------------
# 构造 task
# --------------------------------------------------
tasks <- list()

# Std cycles (8) × modules_uniform (27)
for (i in seq_len(nrow(cycles_std))) {
  for (mod in modules_uniform) {
    tasks[[length(tasks)+1]] <- list(
      filename = paste0(mod, cycles_std$suffix[i]),
      year = cycles_std$pathYear[i], kind = "xpt"
    )
  }
}

# Pre-pandemic P_ modules
pp_modules <- c("DEMO","BMX","BPX","SMQ","ALQ","DIQ","BPQ","MCQ","PAQ","RHQ",
                "PBCD","LUX","BIOPRO","CBC","HSCRP",
                "GHB","GLU","INS","HDL","TCHOL","TRIGLY",
                "DR1TOT","DR2TOT","VID","ALB_CR","COT")
for (mod in pp_modules) {
  tasks[[length(tasks)+1]] <- list(
    filename = paste0("P_", mod),
    year = 2017, kind = "xpt"
  )
}

# PFAS 9 个文件
for (i in seq_len(nrow(pfas_files))) {
  tasks[[length(tasks)+1]] <- list(
    filename = pfas_files$filename[i],
    year = pfas_files$pathYear[i], kind = "xpt"
  )
}

# NCHS LMF for PFAS C-J cycles
mort_base <- "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality"
mort_cycles <- c(2003, 2005, 2007, 2009, 2011, 2013, 2015, 2017)
for (yr in mort_cycles) {
  fn <- paste0("NHANES_", yr, "_", yr+1, "_MORT_2019_PUBLIC.dat")
  tasks[[length(tasks)+1]] <- list(
    filename = fn, year = yr, kind = "mort",
    url = paste0(mort_base, "/", fn)
  )
}

log_line(sprintf("总任务数: %d (xpt=%d, mort=%d)",
                 length(tasks),
                 sum(sapply(tasks, function(t) t$kind == "xpt")),
                 sum(sapply(tasks, function(t) t$kind == "mort"))))

# --------------------------------------------------
# Worker
# --------------------------------------------------
download_task <- function(task, raw_dir, mort_dir) {
  is_real_xpt <- function(path) {
    if (!file.exists(path)) return(FALSE)
    if (file.size(path) < 1024) return(FALSE)
    con <- file(path, "rb"); on.exit(close(con))
    identical(readChar(con, 6, useBytes = TRUE), "HEADER")
  }
  if (task$kind == "xpt") {
    destfile <- file.path(raw_dir, paste0(task$filename, ".xpt"))
    if (is_real_xpt(destfile)) {
      return(list(file = task$filename, status = "已存在"))
    }
    urls <- c(
      sprintf("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/%d/DataFiles/%s.xpt",
              task$year, task$filename),
      sprintf("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/%d/DataFiles/%s.XPT",
              task$year, task$filename)
    )
    for (url in urls) {
      ok <- tryCatch({
        download.file(url, destfile = destfile, mode = "wb",
                      quiet = TRUE, method = "libcurl")
        is_real_xpt(destfile)
      }, error = function(e) FALSE, warning = function(w) FALSE)
      if (isTRUE(ok)) return(list(file = task$filename, status = "ok"))
      if (file.exists(destfile)) file.remove(destfile)
    }
    return(list(file = task$filename, status = "FAIL"))
  } else {
    destfile <- file.path(mort_dir, task$filename)
    if (file.exists(destfile) && file.size(destfile) > 0) {
      return(list(file = task$filename, status = "已存在"))
    }
    ok <- tryCatch({
      download.file(task$url, destfile = destfile, mode = "wb",
                    quiet = TRUE, method = "libcurl")
      file.exists(destfile) && file.size(destfile) > 0
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (isTRUE(ok)) return(list(file = task$filename, status = "ok"))
    return(list(file = task$filename, status = "FAIL"))
  }
}

# --------------------------------------------------
# 并行下载 (8 worker)
# --------------------------------------------------
library(parallel)
n_workers <- min(8L, length(tasks))
log_line(sprintf("启动 %d worker 并行下载 (CDC 限 ~10)...", n_workers))

cl <- makeCluster(n_workers)
on.exit(stopCluster(cl), add = TRUE)

t0 <- Sys.time()
results <- parLapplyLB(cl, tasks, download_task,
                       raw_dir = raw_dir, mort_dir = mort_dir)
elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
stopCluster(cl)

status_tbl <- table(sapply(results, function(r) r$status))
log_line(sprintf("\n并行下载结束 (%.1f 秒):", elapsed))
for (s in names(status_tbl)) {
  log_line(sprintf("  %-10s : %d", s, status_tbl[[s]]))
}

fails <- sapply(results[sapply(results, function(r) r$status == "FAIL")],
                function(r) r$file)
if (length(fails) > 0) {
  log_line("\n失败 (部分模块周期不存在是正常的, 早周期 HSCRP/COT/HEPB_S 等):")
  for (f in fails) log_line(sprintf("  - %s", f))
}

n_xpt  <- length(list.files(raw_dir,  pattern = "\\.xpt$"))
n_mort <- length(list.files(mort_dir, pattern = "\\.dat$"))
log_line(sprintf("\ndata/raw/         %d 个 .xpt", n_xpt))
log_line(sprintf("data/raw/mortality/ %d 个 .dat", n_mort))
log_line("========================================")