# 20e — figures for the CORRECTED paper (PFNA the one independent compound).
# Reads authoritative CSVs (table_xsec_pfas_fib4.csv, table_mutual_pfas_fib4.csv) so
# figures match the manuscript exactly — no in-script refitting.
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
dir.create("output/figures", showWarnings=FALSE, recursive=TRUE)
clean <- function(x){ x<-gsub("_z$","",x); dplyr::recode(x, LBXPFNA="PFNA",LBXPFOS="PFOS",LBXPFOA="PFOA",LBXPFHS="PFHxS",LBXPFDE="PFDA",LBXMPAH="Me-PFOSA-AcOH",.default=x) }
sv <- function(p,f,w=8,h=5.5) ggsave(file.path("output/figures",f),p,width=w,height=h,dpi=300,bg="white")
th <- theme_minimal(base_size=13)+theme(panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=14),axis.text=element_text(color="black"),plot.subtitle=element_text(size=10))

# ---- Fig 1 CONSORT ----
st <- data.frame(y=9:1, n=c(95872,54022,52994,15403,12691,12224,10808,9701,8590),
  lab=c("Pooled NHANES PFAS subsamples, 2005-March 2020","Aged >= 20 years","Not pregnant",
        ">= 1 of six core serum PFAS measured","HBV/HCV-negative","No self-reported liver disease",
        "Complete FIB-4 components","Excl. 2003-04 cycle + 2017-18 standalone de-duplicated",
        "Primary analytic cohort (complete covariates)"))
st$txt <- sprintf("%s\n(n = %s)", st$lab, formatC(st$n,big.mark=",",format="d"))
excl <- c(NA,41850,1028,37591,2712,467,1416,1107,1111)
ar <- data.frame(y=9:2, xl=sprintf("- %s", formatC(excl[2:9],big.mark=",",format="d")))
p1 <- ggplot(st, aes(1,y)) +
  geom_segment(data=st[-9,], aes(x=1,xend=1,y=y-0.62,yend=y-0.40), arrow=arrow(length=unit(0.16,"cm"),type="closed"), color="grey45") +
  geom_text(data=ar, aes(x=1.32, y=y-0.5, label=xl), size=2.9, color="grey35", hjust=0) +
  geom_label(aes(label=txt), fill="grey96", label.size=0.3, size=3.1, lineheight=0.95, label.padding=unit(0.4,"lines")) +
  scale_y_continuous(limits=c(0.4,9.6))+scale_x_continuous(limits=c(0.45,1.75)) +
  labs(title="Figure 1. Participant selection")+theme_void(base_size=12)+theme(plot.title=element_text(face="bold",hjust=0.02,size=13))
sv(p1,"fig1_consort.png",7.5,9)

xs <- read.csv("output/tables/table_xsec_pfas_fib4.csv", stringsAsFactors=FALSE)

# ---- Fig 2 single-PFAS continuous forest (M2) ----
d2 <- xs %>% filter(outcome=="logFIB4", model=="M2") %>%
  mutate(name=clean(pfas), pct=pct_per_SD, lo2=(exp(lo)-1)*100, hi2=(exp(hi)-1)*100,
         sig=ifelse(p<0.05 & est>0,"Positive (P < 0.05)", ifelse(p<0.05,"Inverse (P < 0.05)","Not significant"))) %>% arrange(pct)
d2$name <- factor(d2$name, levels=d2$name)
p2 <- ggplot(d2, aes(pct,name,color=sig)) + geom_vline(xintercept=0,linetype=2,color="grey50") +
  geom_pointrange(aes(xmin=lo2,xmax=hi2),linewidth=0.8,size=0.7) +
  geom_text(aes(label=sprintf("%+.1f%%",pct)),vjust=-0.9,size=3.4,show.legend=FALSE) +
  scale_color_manual(values=c("Positive (P < 0.05)"="#b2182b","Inverse (P < 0.05)"="#2166ac","Not significant"="grey55"),name=NULL) +
  labs(title="Single-PFAS associations with FIB-4 (continuous)",
       subtitle="Per-SD change in FIB-4, fully adjusted (M2); NHANES 2005-March 2020",
       x="Adjusted change in FIB-4 per 1-SD higher PFAS (%)", y=NULL) + th + theme(legend.position="top")
sv(p2,"fig2_pfas_forest.png",8,5.5)

# ---- Fig 3 mutually-adjusted (from table_mutual CSV) ----
mt <- read.csv("output/tables/table_mutual_pfas_fib4.csv", stringsAsFactors=FALSE)
d3 <- mt %>% mutate(name=clean(pfas), pct=pct_per_SD, lo2=(exp(lo)-1)*100, hi2=(exp(hi)-1)*100,
         sig=ifelse(p<0.05 & est>0,"Positive (joint only)", ifelse(p<0.05,"Inverse","Not significant"))) %>% arrange(pct)
d3$name <- factor(d3$name, levels=d3$name)
p3 <- ggplot(d3, aes(pct,name,color=sig)) + geom_vline(xintercept=0,linetype=2,color="grey50") +
  geom_pointrange(aes(xmin=lo2,xmax=hi2),linewidth=0.8,size=0.7) +
  geom_text(aes(label=sprintf("%+.1f%%",pct)),vjust=-0.9,size=3.4,show.legend=FALSE) +
  scale_color_manual(values=c("Positive (joint only)"="#b2182b","Inverse"="#2166ac","Not significant"="grey55"),name=NULL) +
  labs(title="Mutually-adjusted multi-PFAS model",
       subtitle=sprintf("All six PFAS modelled jointly (n=%s); no compound shows a corroborated independent association", formatC(mt$n[1],big.mark=",",format="d")),
       x="Adjusted change in FIB-4 per 1-SD higher PFAS (%), mutually adjusted", y=NULL) + th + theme(legend.position="top")
sv(p3,"fig3_mutual_pfas.png",8,5.5)

# ---- Fig 4 high-positive >=2.67 ORs (M2) ----
d4 <- xs %>% filter(outcome=="FIB4>=2.67", model=="M2") %>%
  mutate(name=clean(pfas), sig=ifelse(p<0.05 & est>1,"Positive (nominal)", ifelse(p<0.05,"Inverse (P < 0.05)","Not significant"))) %>% arrange(est)
d4$name <- factor(d4$name, levels=d4$name)
p4 <- ggplot(d4, aes(est,name,color=sig)) + geom_vline(xintercept=1,linetype=2,color="grey50") +
  geom_pointrange(aes(xmin=lo,xmax=hi),linewidth=0.8,size=0.7) +
  geom_text(aes(label=sprintf("%.2f",est)),vjust=-0.9,size=3.4,show.legend=FALSE) +
  scale_color_manual(values=c("Positive (nominal)"="#b2182b","Inverse (P < 0.05)"="#2166ac","Not significant"="grey55"),name=NULL) + scale_x_log10() +
  labs(title="High-positive advanced fibrosis (FIB-4 >= 2.67)",
       subtitle="Per-SD odds ratios, fully adjusted; no PFAS is positive (PFHxS and Me-PFOSA-AcOH inverse)",
       x="Odds ratio per 1-SD higher PFAS (log scale)", y=NULL) + th + theme(legend.position="top")
sv(p4,"fig4_threshold_or.png",8,5.5)

# ---- Spearman correlations among the six PFAS ----
load("data/processed/nhanes_final.RData")
pf6r <- intersect(c("LBXPFNA_z","LBXPFOS_z","LBXPFOA_z","LBXPFHS_z","LBXPFDE_z","LBXMPAH_z"), names(nhanes_final))
cm <- cor(nhanes_final[,pf6r], use="pairwise.complete.obs", method="spearman")
colnames(cm)<-clean(colnames(cm)); rownames(cm)<-clean(rownames(cm))
cmd <- as.data.frame(as.table(cm)); names(cmd)<-c("x","y","r"); write.csv(cmd, "output/tables/pfas_spearman.csv", row.names=FALSE)
pC <- ggplot(cmd,aes(x,y,fill=r))+geom_tile(color="white")+geom_text(aes(label=sprintf("%.2f",r)),size=3.3)+
  scale_fill_gradient2(low="#2166ac",mid="white",high="#b2182b",midpoint=0,limits=c(-1,1))+
  labs(title="Spearman correlations among serum PFAS",x=NULL,y=NULL)+theme_minimal(base_size=12)+theme(axis.text.x=element_text(angle=45,hjust=1))
sv(pC,"figS_pfas_corr.png",7,6)
cat("Corrected figures written: fig1_consort, fig2_pfas_forest, fig3_mutual_pfas, fig4_threshold_or, figS_pfas_corr\n")