# 20f — NEW graphical abstract (corrected thesis) + Table 1 by PFNA tertile
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(survey) })

# ---------- Graphical abstract (corrective null: crude positive -> adjusted null) ----------
comp  <- c("PFOS","PFOA","PFNA","PFHxS","PFDA","Me-PFOSA-AcOH")
crude <- data.frame(y=6:1, lab=comp)
adj   <- data.frame(y=6:1, lab=comp, grp=c("Null","Null","Null","Null","Null","Inverse"))
pal   <- c("Null"="grey70","Inverse"="#2166ac")
ga <- ggplot() +
  annotate("rect", xmin=0.3, xmax=2.5, ymin=0.4, ymax=6.6, fill="#fde0dd", color="#b2182b") +
  annotate("text", x=1.4, y=7.25, label="Crude: all six PFAS\nassociated with higher FIB-4\n(P < 1e-7)", fontface="bold", size=3.2, lineheight=0.9, color="#b2182b") +
  geom_text(data=crude, aes(x=1.4, y=y, label=lab), size=3.1, fontface="bold") +
  annotate("segment", x=2.7, xend=4.0, y=3.5, yend=3.5, arrow=arrow(length=unit(0.25,"cm"),type="closed"), linewidth=1.1, color="grey30") +
  annotate("text", x=3.35, y=4.5, label="adjust for age, sex,\nrace/ethnicity +\nfull covariates", size=2.9, lineheight=0.9, color="grey25") +
  annotate("rect", xmin=4.2, xmax=6.7, ymin=0.4, ymax=6.6, fill="grey96", color="grey70") +
  annotate("text", x=5.45, y=7.25, label="Adjusted (M2):\nno independent association", fontface="bold", size=3.2, lineheight=0.9, color="grey25") +
  geom_label(data=adj, aes(x=5.45, y=y, label=lab, fill=grp), color="white", fontface="bold", size=3.0, label.padding=unit(0.25,"lines")) +
  scale_fill_manual(values=pal, name=NULL) +
  annotate("label", x=3.5, y=-0.5, label="Crude associations reflect confounding and the age term in FIB-4.\nDirectly measured fibrosis (elastography) shows no PFAS association; mixture and mortality null.", size=2.9, lineheight=1.0, fill="grey92", fontface="italic") +
  scale_x_continuous(limits=c(0,7)) + scale_y_continuous(limits=c(-1.2,7.8)) +
  labs(title="Per- and polyfluoroalkyl substances and liver fibrosis: a compound-resolved NHANES analysis") +
  theme_void(base_size=12) + theme(legend.position="bottom", plot.title=element_text(face="bold", hjust=0.5, size=11)) +
  guides(fill=guide_legend(title=NULL))
ggsave("output/figures/graphical_abstract.png", ga, width=9.5, height=5.4, dpi=300, bg="white")

# ---------- Table 1 by PFNA tertile (survey-weighted) ----------
load("data/processed/nhanes_final.RData"); load("data/processed/nhanes_design.RData")
des <- design_main
des <- update(des, fib4_adv267 = as.integer(exp(fib4_log) >= 2.67),
              pfna_t = cut(LBXPFNA_z, quantile(LBXPFNA_z, c(0,1/3,2/3,1), na.rm=TRUE),
                           labels=c("T1 (low)","T2","T3 (high)"), include.lowest=TRUE))
num_vars <- intersect(c("age","bmi","fbg","hba1c"), names(des$variables))
rows <- list()
rows[["n (unweighted)"]] <- table(des$variables$pfna_t)
for (v in num_vars) {
  m <- svyby(as.formula(paste0("~",v)), ~pfna_t, des, svymean, na.rm=TRUE)
  rows[[paste0(v," (mean)")]] <- sprintf("%.1f", m[[v]])
}
# categorical: % female, % diabetes, FIB-4>=2.67
catpct <- function(var) { t <- svyby(as.formula(paste0("~I(",var,")")), ~pfna_t, des, svymean, na.rm=TRUE); t }
fem <- svyby(~I(RIAGENDR==2), ~pfna_t, des, svymean, na.rm=TRUE)
adv <- svyby(~fib4_adv267, ~pfna_t, des, svymean, na.rm=TRUE)
t1 <- data.frame(
  variable = c("n (unweighted)", paste0(num_vars," (wtd mean)"), "Female (wtd %)", "FIB-4 >= 2.67 (wtd %)"),
  T1 = c(as.integer(table(des$variables$pfna_t)[1]), sapply(num_vars,function(v) sprintf("%.1f", svyby(as.formula(paste0("~",v)),~pfna_t,des,svymean,na.rm=TRUE)[[v]][1])),
         sprintf("%.1f", fem[["I(RIAGENDR == 2)TRUE"]][1]*100), sprintf("%.1f", adv$fib4_adv267[1]*100)),
  T2 = c(as.integer(table(des$variables$pfna_t)[2]), sapply(num_vars,function(v) sprintf("%.1f", svyby(as.formula(paste0("~",v)),~pfna_t,des,svymean,na.rm=TRUE)[[v]][2])),
         sprintf("%.1f", fem[["I(RIAGENDR == 2)TRUE"]][2]*100), sprintf("%.1f", adv$fib4_adv267[2]*100)),
  T3 = c(as.integer(table(des$variables$pfna_t)[3]), sapply(num_vars,function(v) sprintf("%.1f", svyby(as.formula(paste0("~",v)),~pfna_t,des,svymean,na.rm=TRUE)[[v]][3])),
         sprintf("%.1f", fem[["I(RIAGENDR == 2)TRUE"]][3]*100), sprintf("%.1f", adv$fib4_adv267[3]*100)),
  stringsAsFactors=FALSE)
# NOTE: table1_by_pfna.csv is the authoritative full Table 1 written by _table1_full.R (run afterwards).
# 20f produces only the graphical abstract; the short t1 here is for the console preview only.
cat("graphical_abstract.png written; Table 1 by PFNA tertile:\n"); print(t1)