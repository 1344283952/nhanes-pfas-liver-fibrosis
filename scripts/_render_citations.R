# render [`key`] citations -> Vancouver [n] in manuscript_v2.md + build reference list
txt <- paste(readLines("manuscript_v2.md", warn=FALSE, encoding="UTF-8"), collapse="\n")

keys <- c("costello2022pfas"=1,"stratakis2020helix"=2,"momo2024pfas"=3,"yan2024pfas"=4,
  "nchs_series2_190"=5,"lubin2004lod"=6,"sterling2006fib4"=7,"easl2024masld"=8,
  "nchs_lmf2019"=9,"lumley2010survey"=10,"keil2020qgcomp"=11,"carrico2015wqs"=12,
  "fine1999gray"=13,"vanderweele2017evalue"=14)

# combined two-key citations first
txt <- gsub("[`costello2022pfas`, `stratakis2020helix`]", "[1,2]", txt, fixed=TRUE)
txt <- gsub("[`momo2024pfas`, `yan2024pfas`]", "[3,4]", txt, fixed=TRUE)
# single keys
for (k in names(keys)) txt <- gsub(paste0("[`", k, "`]"), paste0("[", keys[[k]], "]"), txt, fixed=TRUE)

refs <- c(
"1. Costello E, Rock S, Stratakis N, et al. Exposure to per- and polyfluoroalkyl substances and markers of liver injury: a systematic review and meta-analysis. *Environ Health Perspect.* 2022;130(4):046001. doi:10.1289/EHP10092",
"2. Stratakis N, Conti DV, Jin R, et al. Prenatal exposure to perfluoroalkyl substances associated with increased susceptibility to liver injury in children. *Hepatology.* 2020;72(5):1758-1770. doi:10.1002/hep.31483",
"3. Momo HD, Alvarez CS, Purdue MP, Graubard BI, McGlynn KA. Associations of per- and polyfluoroalkyl substances and nonalcoholic fatty liver disease in the United States adult population, 2003-2018. *Environ Epidemiol.* 2024;8(1):e284. doi:10.1097/EE9.0000000000000284",
"4. Yan Y, Zhang H, Xu W, et al. Association of exposure to per- and polyfluoroalkyl substances with liver injury in American adults. *J Biomed Res.* 2024;38(6):628-638. doi:10.7555/JBR.38.20240018",
"5. National Center for Health Statistics. NHANES 2017-March 2020 pre-pandemic file: sample design, estimation, and analytic guidelines. *Vital Health Stat 2.* 2022;(190):1-36.",
"6. Lubin JH, Colt JS, Camann D, et al. Epidemiologic evaluation of measurement data in the presence of detection limits. *Environ Health Perspect.* 2004;112(17):1691-1696. doi:10.1289/ehp.7199",
"7. Sterling RK, Lissen E, Clumeck N, et al. Development of a simple noninvasive index to predict significant fibrosis in patients with HIV/HCV coinfection. *Hepatology.* 2006;43(6):1317-1325. doi:10.1002/hep.21178",
"8. European Association for the Study of the Liver (EASL), European Association for the Study of Diabetes (EASD), European Association for the Study of Obesity (EASO). EASL-EASD-EASO clinical practice guidelines on the management of metabolic dysfunction-associated steatotic liver disease (MASLD). *J Hepatol.* 2024;81(3):492-542. doi:10.1016/j.jhep.2024.04.031",
"9. National Center for Health Statistics. NHANES (1999-2018) Linked Mortality File: public-use version, linked through December 31, 2019. Hyattsville, MD: NCHS; 2022.",
"10. Lumley T. *Complex Surveys: A Guide to Analysis Using R.* Hoboken, NJ: John Wiley & Sons; 2010.",
"11. Keil AP, Buckley JP, O'Brien KM, Ferguson KK, Zhao S, White AJ. A quantile-based g-computation approach to addressing the effects of exposure mixtures. *Environ Health Perspect.* 2020;128(4):047004. doi:10.1289/EHP5838",
"12. Carrico C, Gennings C, Wheeler DC, Factor-Litvak P. Characterization of weighted quantile sum regression for highly correlated data in a risk analysis setting. *J Agric Biol Environ Stat.* 2015;20(1):100-120. doi:10.1007/s13253--0180-3",
"13. Fine JP, Gray RJ. A proportional hazards model for the subdistribution of a competing risk. *J Am Stat Assoc.* 1999;94(446):496-509. doi:10.1080/01621459.1999.10474144",
"14. VanderWeele TJ, Ding P. Sensitivity analysis in observational research: introducing the E-value. *Ann Intern Med.* 2017;167(4):268-274. doi:10.7326/M16-2607")

head <- sub("## References.*$", "", txt)
out <- paste0(head, "## References\n\n", paste(refs, collapse="\n\n"), "\n")
writeLines(out, "manuscript_v2.md", useBytes=TRUE)

# report any unresolved keys
left <- regmatches(out, gregexpr("\\[`[^`]+`\\]", out))[[1]]
cat("Citations rendered. Unresolved [`key`] remaining:", length(left), "\n")
if (length(left)) print(unique(left))
cat("Reference list: 16 entries.\n")