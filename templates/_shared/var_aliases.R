
suppressPackageStartupMessages({
  library(dplyr)
})

if (!exists("coalesce_cols", mode = "function")) {
  coalesce_cols <- function(df, cols) {
    cols <- intersect(cols, names(df))
    if (length(cols) == 0) return(rep(NA_real_, nrow(df)))
    if (length(cols) == 1) return(df[[cols[1]]])
    Reduce(function(a, b) ifelse(is.na(a), b, a),
           lapply(cols, function(c) df[[c]]))
  }
}


NHANES_VAR_ALIASES <- list(

  hdl_mgdl     = c("LBDHDL", "LBXHDD", "LBDHDD"),
  tchol_mgdl   = c("LBXTC"),
  tg_mgdl      = c("LBXTR"),
  ldl_mgdl     = c("LBDLDL"),

  ast_unl      = c("LBXSASSI", "LBXSAST"),
  alt_unl      = c("LBXSATSI", "LBXSALT"),
  alp_unl      = c("LBXSAPSI", "LBXSAPSU"),
  ggt_unl      = c("LBXSGTSI"),
  total_bili_mgdl = c("LBXSTB"),
  alb_gl       = c("LBDSALSI"),
  alb_gdl      = c("LBXSAL"),

  creatinine_serum = c("LBXSCR"),
  bun_mgdl     = c("LBXSBU"),
  uric_acid    = c("LBXSUA"),

  hba1c_pct    = c("LBXGH"),
  fpg_mgdl     = c("LBXGLU"),
  insulin_uum  = c("LBXIN"),

  platelet_x10e9 = c("LBXPLTSI"),
  wbc_x10e9      = c("LBXWBCSI"),
  rbc_x10e6_ul   = c("LBXRBCSI"),
  hgb_gdl        = c("LBXHGB"),
  neutrophil_abs = c("LBDNENO"),
  lymphocyte_abs = c("LBDLYMNO"),
  monocyte_abs   = c("LBDMONO"),
  eosinophil_abs = c("LBDEONO"),
  basophil_abs   = c("LBDBANO"),

  lead_ugdl    = c("LBXBPB"),
  cadmium_ugl  = c("LBXBCD"),
  hg_total_ugl = c("LBXTHG"),
  se_whole_blood = c("LBXBSE"),
  manganese_ugl  = c("LBXBMN"),

  hg_inorganic_ugl = c("LBXIHG"),
  hg_methyl_ugl    = c("LBXBGM"),

  uri_arsenic_total = c("LBXTUA"),
  uri_cadmium       = c("URXUCD"),
  uri_lead          = c("URXUPB"),
  uri_mercury       = c("URXUHG"),
  uri_creatinine    = c("URXUCR"),

  se_dietary       = c("DR1TSELE"),
  se_dietary_d2    = c("DR2TSELE"),
  zn_dietary       = c("DR1TZINC"),
  zn_dietary_d2    = c("DR2TZINC"),
  cu_dietary       = c("DR1TCOPP"),
  cu_dietary_d2    = c("DR2TCOPP"),
  kcal_dietary     = c("DR1TKCAL"),
  protein_g        = c("DR1TPROT"),
  carb_g           = c("DR1TCARB"),
  fat_g            = c("DR1TTFAT"),
  fiber_g          = c("DR1TFIBE"),

  cap_median       = c("LUXCAPM"),
  cap_iqr          = c("LUXCAPIQR"),
  lsm_median       = c("LUXSMED"),
  lsm_iqr_med      = c("LUXSIQRM"),
  lux_n_attempt    = c("LUXSREPC"),

  sbp_1 = c("BPXSY1"), sbp_2 = c("BPXSY2"), sbp_3 = c("BPXSY3"), sbp_4 = c("BPXSY4"),
  dbp_1 = c("BPXDI1"), dbp_2 = c("BPXDI2"), dbp_3 = c("BPXDI3"), dbp_4 = c("BPXDI4"),

  bmi              = c("BMXBMI"),
  waist_cm         = c("BMXWAIST"),
  height_cm        = c("BMXHT"),
  weight_kg        = c("BMXWT"),
  hip_cm           = c("BMXHIP"),

  drink_lifetime   = c("ALQ101", "ALQ111"),
  drink_freq       = c("ALQ121", "ALQ120Q"),
  drink_amount     = c("ALQ130", "ALQ142"),

  smoke_ever       = c("SMQ020"),
  smoke_now        = c("SMQ040"),
  smoke_age_start  = c("SMD030"),
  smoke_cigs_per_day = c("SMD650", "SMD641"),

  ever_pregnant    = c("RHQ160", "RHD143"),
  gdm_history      = c("RHQ162"),
  parity           = c("RHQ171", "RHD180"),
  age_first_period = c("RHQ"),

  sleep_hours      = c("SLDH", "SLD"),
  sleep_disorder   = c("SLQ050"),
  sleep_troubled   = c("SLQ120"),

  diabetes_self    = c("DIQ"),
  diabetes_insulin = c("DIQ050"),
  diabetes_oralmed = c("DIQ070"),
  htn_self         = c("BPQ020"),
  htn_med          = c("BPQ050A"),
  lipid_med        = c("BPQ090D", "BPQ100D"),

  liver_disease_self = c("MCQ160L"),
  liver_cancer_self  = c("MCQ240I", "MCQ240L"),
  cancer_ever        = c("MCQ220"),

  family_hx_heart    = c("MCQ300A"),
  family_hx_diabetes = c("MCQ300C"),

  vit_d_serum      = c("LBXVIDMS", "LBXVDMS", "LBXVD2MS", "LBXVD3MS"),

  hbsag = c("LBXHBS"),
  hcv   = c("LBXHCR")
)


apply_var_aliases <- function(df, alias_list = NHANES_VAR_ALIASES, skip_existing = TRUE) {
  for (canonical in names(alias_list)) {
    if (skip_existing && canonical %in% names(df)) next
    candidates <- alias_list[[canonical]]
    avail <- intersect(candidates, names(df))
    if (length(avail) == 0) next  # no candidate columns in df; skip
    df[[canonical]] <- coalesce_cols(df, candidates)
  }
  df
}


report_alias_coverage <- function(df, alias_list = NHANES_VAR_ALIASES) {
  out <- lapply(names(alias_list), function(canonical) {
    candidates <- alias_list[[canonical]]
    avail <- intersect(candidates, names(df))
    if (canonical %in% names(df)) {
      n_nn <- sum(!is.na(df[[canonical]]))
      pct  <- 100 * n_nn / nrow(df)
    } else {
      n_nn <- 0L; pct <- 0
    }
    data.frame(
      canonical = canonical,
      candidates_present = paste(avail, collapse = ","),
      n_nonNA = n_nn,
      pct_nonNA = round(pct, 1),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}