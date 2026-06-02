# templates/_consistency_check.R
# 数字一致性扫描 v1（投稿前 QA Step B 内）
#
# 目的：抽出 manuscript 系列 .md 文件里所有 HR / OR / aHR / RR / β / CI / P
# 值，导出 CSV 留 PI / 学生人工复核（v2 才做自动 cross-table-vs-manuscript 比对）。
#
# 用法（必须 cd 到 projects/<NNN>/）：
#   Rscript ../../templates/_consistency_check.R
#
# 输出：
#   output/_consistency_report.csv
#
# 退出码：
#   0 = 抽取成功（不代表数字一致，只是抽出了）
#   1 = 找不到任何 manuscript 文件 / 抽到 0 数字
#
# 不做什么（留 v2 / PI 手工）：
#   - 不自动比对 manuscript 数字 vs output/tables/*.xlsx 的真值
#   - 不判断哪些数字是 "应该一致" 的（如 abstract HR vs Results HR）

suppressPackageStartupMessages({
  library(stringr)
})

manuscript_files <- c(
  "manuscript.md", "manuscript.docx",
  "结果说明.md", "方法学段.md",
  "cover_letter.md"
)
report_path <- "output/_consistency_report.csv"

# 抽数字模式：metric = value 格式
patterns <- list(
  hr_or  = "\\b(aHR|HR|aOR|OR|RR|HR_adj)\\b\\s*[=＝]\\s*(-?\\d+\\.\\d{2,3})",
  ci     = "(95% ?CI)[: ]*(-?\\d+\\.\\d{2,3})\\s*[-–]\\s*(-?\\d+\\.\\d{2,3})",
  pval   = "\\b([Pp](?:-trend|-interaction|_int)?)\\s*[=＝<>]\\s*(\\d?\\.\\d{2,5})",
  beta   = "\\b(beta|β)\\s*[=＝]\\s*(-?\\d+\\.\\d{2,3})",
  evalue = "\\b(E-value[s]?)\\s*[=＝]\\s*(-?\\d+\\.\\d{1,2})"
)

extract_one_file <- function(path) {
  if (!file.exists(path)) return(NULL)
  txt <- tryCatch(
    paste(readLines(path, encoding = "UTF-8", warn = FALSE), collapse = "\n"),
    error = function(e) NULL
  )
  if (is.null(txt) || nchar(txt, type = "bytes") == 0) return(NULL)
  rows <- list()
  for (type in names(patterns)) {
    mat <- str_match_all(txt, patterns[[type]])[[1]]
    if (nrow(mat) == 0) next
    for (i in seq_len(nrow(mat))) {
      rows[[length(rows) + 1L]] <- data.frame(
        source = basename(path),
        type   = type,
        metric = if (ncol(mat) >= 2) mat[i, 2] else NA_character_,
        value  = if (ncol(mat) >= 3) mat[i, 3] else mat[i, 2],
        value2 = if (ncol(mat) >= 4) mat[i, 4] else NA_character_,
        full   = mat[i, 1],
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

main <- function() {
  found <- list()
  for (f in manuscript_files) {
    res <- extract_one_file(f)
    if (!is.null(res) && nrow(res) > 0) found[[f]] <- res
  }

  if (length(found) == 0) {
    message("[warn] no manuscript files matched any of: ",
            paste(manuscript_files, collapse = ", "))
    message("[warn] cwd = ", getwd())
    return(invisible(1L))
  }

  out <- do.call(rbind, found)
  dir.create(dirname(report_path), showWarnings = FALSE, recursive = TRUE)
  write.csv(out, report_path, row.names = FALSE, fileEncoding = "UTF-8")

  message(sprintf("[ok] %d numeric claims extracted from %d files -> %s",
                  nrow(out), length(found), report_path))
  message("[next] open CSV in Excel; sort by metric+value;")
  message("       same metric appearing with different values across files = inconsistency")
  message("       cross-check against output/tables/*.xlsx for the underlying truth")

  invisible(0L)
}

if (sys.nframe() == 0L) {
  status <- main()
  if (!is.null(status)) quit(save = "no", status = status)
}