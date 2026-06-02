# nhanes-pfas-liver-fibrosis

> **Pre-publication notice — please read.** This repository accompanies a research manuscript that is **currently under peer review and not yet published**. The study, its results, and its figures have not been peer-reviewed and remain the authors' unpublished work. The code is shared under the MIT License for transparency and reproducibility; however, **please do not reproduce, reuse, or build upon the results, figures, or analysis pipeline in any other publication, preprint, thesis, or grant application until this article is published.** For any use, please cite the published article (once available) and contact the corresponding author, Ling Yu (yulingyxb@jlu.edu.cn). Thank you for respecting the authors' priority.

Reproducibility package for a compound-resolved, multi-index NHANES analysis of per- and polyfluoroalkyl substances (PFAS) and liver fibrosis (FIB-4, the age-free APRI, and directly measured transient elastography), with all-cause and cardiometabolic mortality as secondary outcomes.

## Repository layout

```
scripts/             analysis pipeline + helpers (run_all.R orchestrates the full run)
templates/_shared/   shared functions sourced by the pipeline (FIB-4, survey weights, covariates)
data/processed/      processed analytic datasets (*.RData) -- included
data/raw/            empty; NHANES XPT files are downloaded here by scripts/01_download_data.R
output/tables/       result tables backing the manuscript
output/figures/      manuscript figures (300 DPI PNG)
```

## Full pipeline

```r
Rscript scripts/00_install_packages.R   # one-time: install R package dependencies
Rscript scripts/run_all.R               # download (~700 MB) -> clean -> survey design -> analysis -> figures
```

## Inspect results without re-downloading

The processed analytic datasets (`data/processed/*.RData`) and every committed file under `output/` are the outputs of the pipeline, so the reported tables and figures can be inspected directly. To re-run only the analysis on the included data (no ~700 MB download), run the analysis scripts after `00_install_packages.R` in the order listed in `scripts/run_all.R` (from `06d_reanalysis_cleanM2.R` onward).

## Notes

- Raw NHANES XPT files are public-domain and not vendored; re-download with `scripts/01_download_data.R`.
- All models use the PFAS-subsample serum weights and the NHANES complex-survey design.
- The mixture analysis uses quantile g-computation (`qgcomp`) and weighted quantile sum regression.
- Before any analysis, `templates/_shared/_integrity_preflight.R` asserts the AST/ALT mapping and FIB-4 formula direction (guards against a silently inverted fibrosis index).

## License

MIT (see LICENSE).
