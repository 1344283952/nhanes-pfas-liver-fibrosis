
.find_shared <- function(fname) {
  cands <- c(
    file.path("templates", "_shared", fname),          # 从仓库根
    file.path("..", "..", "templates", "_shared", fname), # 从 projects/<NNN>/
    file.path(getwd(), fname),                          # 与本脚本同目录
    fname
  )
  hit <- cands[file.exists(cands)]
  if (length(hit) == 0)
    stop("_integrity_preflight: 找不到 ", fname,
         " —— 尝试过: ", paste(cands, collapse = " | "))
  hit[1]
}

source(.find_shared("var_aliases.R"))
source(.find_shared("fib4.R"))

fail <- character(0)

ast1 <- NHANES_VAR_ALIASES$ast_unl[1]
alt1 <- NHANES_VAR_ALIASES$alt_unl[1]
if (!identical(ast1, "LBXSASSI"))
  fail <- c(fail, sprintf("AST 映射错: ast_unl[1]='%s' 应为 'LBXSASSI' (CDC: LBXSASSI=AST)", ast1))
if (!identical(alt1, "LBXSATSI"))
  fail <- c(fail, sprintf("ALT 映射错: alt_unl[1]='%s' 应为 'LBXSATSI' (CDC: LBXSATSI=ALT)", alt1))
if ("LBXSATSI" %in% NHANES_VAR_ALIASES$ast_unl)
  fail <- c(fail, "LBXSATSI(实为 ALT)出现在 ast_unl 候选里 —— 互换风险")
if ("LBXSASSI" %in% NHANES_VAR_ALIASES$alt_unl)
  fail <- c(fail, "LBXSASSI(实为 AST)出现在 alt_unl 候选里 —— 互换风险")

got  <- calc_fib4(age = 60, ast = 80, alt = 20, plt = 150)
want <- (60 * 80) / (150 * sqrt(20))   # 正确公式硬编码 = 7.1554; 互换会得 0.894
if (!isTRUE(abs(got - want) < 1e-6))
  fail <- c(fail, sprintf(
    "FIB-4 公式方向错: calc_fib4(60,80,20,150)=%.4f 应=%.4f (AST 须在分子, ALT 在分母 sqrt)",
    got, want))

if (length(fail) > 0) {
  cat("\n!!!!!!!!!! 完整性预检失败 (INTEGRITY PRE-FLIGHT FAILED) !!!!!!!!!!\n")
  for (f in fail) cat("  [FAIL] ", f, "\n", sep = "")
  cat("\n修复 templates/_shared/var_aliases.R 与 fib4.R 后重跑。\n")
  cat("背景见 memory: feedback_fib4_ast_alt_swap_template_bug。\n")
  cat("⚠️ 任何在修复前跑过 03_clean_data.R 的项目 (含 this project) 都须重跑复核。\n\n")
  stop("integrity pre-flight failed —— 已中止, 以防产出 AST/ALT 互换的假结果")
}

cat(sprintf(
  "[OK] 完整性预检通过: AST=LBXSASSI / ALT=LBXSATSI 映射正确; FIB-4 公式方向正确 (smoke=%.3f)\n",
  got))