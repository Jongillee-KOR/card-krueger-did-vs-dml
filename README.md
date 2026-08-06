# Card & Krueger Replication: Classic DiD vs. Double Machine Learning

**Authors:** Jonggil Lee and Jongju Lee

Using the original Card & Krueger (1994) survey of NJ and PA fast-food restaurants around New Jersey's 1992 minimum wage increase, this project asks whether the classic difference-in-differences estimate of the employment effect holds up once we swap out ordinary least squares for a more flexible, machine-learning-based estimator — and then goes a step further, stress-testing an unexpected pattern in the standard errors that showed up once the two methods were put on genuinely equal footing.

## Motivation

Card & Krueger's NJ-PA fast-food study has been replicated more times than almost any other paper in labor economics, so reproducing the original numbers isn't really the interesting part anymore. What's more useful to check is whether the result holds up when the estimation method changes, not just the data.

Card & Krueger's original finding, that raising the minimum wage didn't cost jobs and if anything nudged employment up slightly, was controversial partly because it came from a fairly simple linear regression. If that result depended heavily on functional form (a linear model, specific controls), a more flexible estimator might tell a different story.

That's what we wanted to check here. We estimate the same effect two ways: a classic OLS difference-in-differences (matching the original paper's approach), and a double/debiased machine learning (DML) estimator that lets Lasso and XGBoost flexibly model the relationship between the covariates and both the treatment (the wage gap) and the outcome (employment change), then partials that out before estimating the effect. If both approaches land in the same place, that's real evidence the original result isn't an artifact of the linear specification.

The two approaches, side by side:

```
ΔFTE = α + β(gap) + ε                    (classic DiD, OLS)

Ỹ = Y − E[Y|X],  D̃ = D − E[D|X]          (Robinson (1988) partialling-out)
Ỹ = θ(D̃) + ε                             (DML: θ estimated via cross-fitted ML for E[Y|X], E[D|X])
```

This started as a basic classic-DiD replication of Card & Krueger's NJ vs. PA fast-food employment study. Everything after that — the robustness check, the DML estimation, the final comparison, and the deeper investigation into why the two methods' standard errors behaved unexpectedly — we built together.

## Key Results

**Classic DiD vs. DML: gap coefficient (primary comparison, unified sample N=370, HC1-robust SE)**

| Method | Estimate | SE | 95% CI |
|---|---|---|---|
| Classic DiD (OLS, unified sample, HC1-robust) | 14.10 | 7.13 | (0.12, 28.08) |
| DML (Lasso) | 11.64 | 6.98 | (-2.04, 25.33) |
| DML (XGBoost) | 11.68 | 6.85 | (-1.74, 25.10) |

**NJ vs. PA employment change (for reference, the original Card-Krueger comparison)**

| | PA | NJ | NJ − PA |
|---|---|---|---|
| FTE before | 23.33 (SE 1.35) | 20.44 (SE 0.51) | -2.89 (SE 1.44) |
| FTE after | 21.17 (SE 0.94) | 21.03 (SE 0.52) | -0.14 (SE 1.08) |
| Change | -2.17 | +0.59 | **+2.75 (SE 1.80)** |

**Interpretation**

The first table is this project's own contribution, and it's the "unified" version of the comparison: Classic DiD and both DML learners are now estimated on the exact same sample (N=370, missing chain/ownership values encoded as their own category rather than dropped) and with the same SE estimator (HC1-robust), so this is a genuine apples-to-apples comparison in a way the original three-way comparison wasn't (see Limitations for how that gap used to work). All three estimates land in a fairly tight band (11.64 to 14.10), and swapping a linear model for two very different flexible ML methods still doesn't move the point estimate very far. That's a meaningful robustness check: if the OLS result had been an artifact of assuming a linear relationship between the covariates and the outcome, we'd have expected Lasso or XGBoost to pull the estimate somewhere else. They didn't.

One thing that did change once the sample and SE estimator were unified: Classic DiD's 95% CI (0.12, 28.08) is the only one of the three that excludes zero, if only barely, while both DML CIs still cross it. That's a real but fragile difference — a lower bound of 0.12 is not something to lean on hard — and it's worth reading alongside the deeper standard-error investigation below, since the way we measure Classic DiD's SE (asymptotic HC1 vs. bootstrap) turns out to matter more than we initially expected.

The second table is just the baseline the comparison above rests on, and it reproduces Card & Krueger's headline result almost exactly: PA lost jobs (-2.17 average FTE) after the law took effect, while NJ actually gained a bit (+0.59), for a NJ-PA gap of +2.75. It's the well-known finding, included here mainly to confirm the replication is sound before trusting the method comparison built on top of it.

## Digging Into the SE Reversal (Scripts 07-09)

When we moved Classic DiD onto the same unified sample and HC1-robust SE as DML, something unexpected showed up: Classic DiD's standard error (7.13) came out *larger* than both DML-Lasso's (6.98) and DML-XGBoost's (6.85). That's backwards from the textbook prediction — in a low-dimensional setting like this one, OLS is supposed to be the more efficient (lower-variance) estimator, with DML's flexibility typically costing some efficiency in exchange for robustness to functional form. We built three follow-up scripts (07-09) to figure out whether this "reversal" was a real methodological finding or an artifact of something narrower in this particular sample.

**Working hypothesis: sparse cells.** Our first guess was that the 8-cell chain × ownership interaction implied by OLS's dummy encoding — with the smallest cell (Chain_4 × Co_Owned) holding only 10 observations — was disadvantaging OLS relative to DML's regularized first stage, which handles thin cells more gracefully.

**07.1 — Coarser categorical encoding.** We collapsed chain from 5 levels down to 2 (the most common chain vs. everything else), which folds the sparse cell into a much larger "Other" group. If the sparse-cell story were the whole explanation, Classic DiD's SE should have dropped below the DML SEs. It barely moved at all (7.13 → 7.13, a change in the fourth decimal place) and never dropped below either DML SE (13.39 for Lasso, 6.99 for XGBoost) — at best, partial support for the hypothesis.

**07.2 — Bootstrap SE.** Instead of asymptotic HC1 SEs, we resampled the data directly (1,000 reps for Classic DiD, 500 for each DML learner, refitting the full cross-fitting procedure each time for DML). Under bootstrap, the reversal disappears entirely: Classic DiD's SE (7.31) becomes the *smallest* of the three, with DML-Lasso at 7.33 and DML-XGBoost at 7.46 — the textbook ordering is restored. This points to the asymptotic HC1 SE being the less reliable estimator here, not to Classic DiD genuinely being less efficient.

**07.3 — Seed sensitivity.** Since DML's cross-fitting depends on a random 5-fold split, we re-ran it across 30 different seeds. DML-Lasso's SE never exceeded Classic DiD's fixed HC1 SE (0 of 30 seeds; range 6.82-7.01), and DML-XGBoost's SE exceeded it in only 2 of 30 seeds (range 6.69-7.18). So under the asymptotic HC1 framing, the reversal is not just a lucky fold split — it's fairly consistent.

**08 — Sparse-cell exclusion.** As a more direct test than the coarsening above, we dropped just the 10 sparse-cell rows and re-ran everything on the remaining N=360. If that one cell were driving the reversal, dropping it should have flipped the pattern. It didn't: the HC1 reversal (Classic DiD SE largest) was still present in the excluded sample, and the bootstrap non-reversal (Classic DiD SE smallest) was also still present. The reversal pattern is unchanged either way, so the single sparse cell isn't the driver — the choice of SE estimator matters far more than this one cell.

**09 — Monte Carlo simulation with a known true effect.** Scripts 07-08 only test the observed data — a single real-world draw. To rule out the whole pattern being a fluke of this one sample, we built a synthetic DGP that resamples real rows (preserving the actual chain × ownership joint distribution, sparse cell included) but assigns a known true effect (β = 12) plus fresh calibrated noise, then ran 500 replications. The results: all three methods have small bias (0.15-0.36, tiny relative to a standard deviation of roughly 7), and no method uniformly dominates on variance — empirical SD is 6.89 for Classic DiD, 6.87 for DML-Lasso, and 7.19 for DML-XGBoost. More importantly, all three methods' *own reported* HC1 SE understates their true (empirical) sampling variability by a similar margin (Classic DiD: 6.80 reported vs. 6.89 actual; Lasso: 6.70 vs. 6.87; XGBoost: 6.85 vs. 7.19) — yet 95% coverage still lands close to nominal for all three (94.4%-94.6%). The miscalibration is real but modest, and it isn't specific to Classic DiD.

**Bottom line.** Put together, 07-09 suggest the SE reversal we saw when unifying the sample is less a story about DML being more efficient than OLS, and more a story about asymptotic HC1 standard errors being a somewhat shaky estimator at this sample size and this covariate structure — for all three methods, not something unique to Classic DiD. Switch to bootstrap SE and the textbook ordering (OLS SE ≤ DML SE) reappears; the underlying point estimates and relative efficiency, confirmed by the Monte Carlo check, are genuinely close across all three methods.

## Limitations & Future Direction

Sample size and missing data: the primary complete-case estimation sample only keeps stores with complete FTE and wage data in both waves (351 of the 410 surveyed), and the robustness-check sample that imputes missing chain/ownership values as their own "Missing" category brings that up to 370. Neither is the full survey, and if the missingness isn't random, both samples could be biased in ways we haven't tested for.

**The Classic DiD vs. DML comparison is now genuinely apples-to-apples.** Earlier in this project, Classic DiD (model5) was estimated on the complete-case sample (N=351, classical SE) while both DML models ran on the imputed sample (N=370, HC1-robust SE) — so the comparison was changing both the sample and the SE estimator between methods, not just the method itself. That complete-case, classical-SE estimate (12.78, SE 6.67) is still saved for reference, but the primary comparison used throughout this README now refits Classic DiD on the identical N=370 sample with the identical HC1-robust SE that DML uses, closing that gap.

**Asymptotic HC1 SE is not fully trustworthy at this sample size.** The 07-09 investigation found that asymptotic HC1 SEs are somewhat unreliable for all three methods in this small, sparse-cell sample — not uniquely for Classic DiD. Bootstrap SEs give a materially different (and, per the Monte Carlo check, better-calibrated) picture, and even change which method looks more "efficient." Any single point estimate/SE pair from this pipeline should be read with that caveat in mind.

**Wide, and occasionally fragile, confidence intervals.** All three methods produce standard errors in the 6.8-7.1 range under HC1, wide enough that even a fairly large point estimate (11.6-14.1) barely clears conventional significance thresholds, if it clears them at all. Classic DiD's CI under the primary specification (0.12, 28.08) technically excludes zero, but only by a hair, and the bootstrap SE for the same coefficient would widen that interval enough to cross zero again — so this shouldn't be read as strong evidence against a null effect.

DML's covariate set here is basically the same as the OLS model's (chain, ownership, region), just modeled more flexibly. That's enough to show the linear DiD result isn't an artifact of functional form, but it doesn't test whether a genuinely richer set of controls (initial wage levels, initial FTE, survey timing) would move the estimate further.

**Next steps.** A cluster or wild bootstrap tailored specifically to the cross-fitting procedure (rather than the simple case-resampling bootstrap used in 07-08) would sharpen the DML SE estimates further. It would also be worth testing a genuinely richer covariate set (beyond chain/ownership/region) to see whether DML's flexibility starts to pay off once there's more nonlinear structure for it to find, and repeating the Monte Carlo check (09) with a couple of different `true_beta` values and DGP shapes to see how sensitive the bias/coverage findings are to the specific calibration chosen here.

## Data

- **Original survey data**: Card & Krueger's public-use NJ/PA fast-food restaurant survey (`public.dat`, `codebook`), covering wave 1 (February/March 1992, before the NJ minimum wage increase) and wave 2 (November/December 1992, after).
- **FTE employment**: constructed as full-time employees + number of managers + 0.5 × part-time employees, separately for each wave.
- **Treatment (`gap`)**: the percentage wage increase a NJ store needed to reach the new $5.05 minimum; zero for PA stores and for NJ stores already paying above $5.05.

All series come from the same original dataset; nothing external is merged in.

## Methodology

| Step | Script | What it does |
|---|---|---|
| 0 | `00_setup.R` | Loads/installs required packages, creates output directories |
| 1 | `01_data_cleaning.R` | Parses the codebook, builds FTE and wage-gap variables, constructs the primary complete-case estimation sample |
| 2 | `02_summary_tables.R` | Reproduces the original paper's summary tables: mean FTE before/after by state, and the NJ-PA comparison |
| 3 | `03_classic_did_regression.R` | Primary analysis: OLS DiD with the NJ dummy and continuous gap specifications, nested F-tests for chain/ownership/region controls; also refits model5 with HC1-robust SE and on the DML-unified sample (N=370) so the method comparison is apples-to-apples |
| 4 | `04_robustness check.R` | Re-runs the same regressions on a larger sample that imputes missing values as zero (and flags what was originally missing) instead of dropping them |
| 5 | `05_dml.R` | Cross-fitted double/debiased ML (Lasso and XGBoost nuisance models, 5-fold cross-fitting) on the larger imputed sample |
| 6 | `06_comparison.R` | Builds the final comparison table and heatmap: Classic DiD vs. DML (Lasso) vs. DML (XGBoost) |
| 7 | `07_Robustness Check.R` | Tests the sparse-cell hypothesis for the SE reversal: coarser chain encoding, bootstrap SE, and fold-seed sensitivity |
| 8 | `08_SparseCell Exclusion Check.R` | Drops the 10-row sparse cell (Chain_4 × Co_Owned) and re-runs the HC1 and bootstrap comparison on N=360, to measure that cell's actual influence |
| 9 | `09_Simulation Check.R` | Monte Carlo simulation (500 reps) with a known true effect on a synthetic DGP that preserves the real covariate structure, to validate the 07-08 findings against ground truth |

## Figures

![Comparison of Minimum Wage Impact: Classic DiD vs. DML](output/figures/classic_did_vs_dml.png)

![Model Comparison Heatmap](output/figures/comparison_heatmap.png)

- `classic_did_vs_dml.png` — point estimate and 95% CI, all three methods side by side (shown above)
- `comparison_heatmap.png` — estimate and SE across all three methods (shown above)
- `dml_lasso_residuals.png` — orthogonalized treatment vs. outcome residuals, Lasso nuisance models
- `dml_xgboost_residuals.png` — orthogonalized treatment vs. outcome residuals, XGBoost nuisance models
- `coarser_categorical_comparison.png` — 5-level vs. 2-level chain encoding, estimate & CI for all three methods (07.1)
- `bootstrap_se_comparison.png` — bootstrap distribution of the gap coefficient, all three methods (07.2)
- `seed_sensitivity_boxplot.png` — DML SE across 30 fold-assignment seeds vs. Classic DiD's fixed SE (07.3)
- `sparse_cell_exclusion_comparison.png` — HC1 and bootstrap SE, full sample vs. sparse-cell-excluded sample (08)
- `simulation_check_boxplot.png` — Monte Carlo distribution of estimates vs. the known true effect (09)

## Repository Structure

```
Classic DiD vs DML/
├── 00_setup.R                        # loads/installs required packages, creates output folders
├── 01_data_cleaning.R                # parses codebook, builds FTE/gap variables, complete-case sample
├── 02_summary_tables.R               # reproduces original paper's Tables 1-3
├── 03_classic_did_regression.R       # primary OLS DiD analysis + nested F-tests + unified-sample HC1-robust refit
├── 04_robustness check.R             # imputed-sample robustness check
├── 05_dml.R                          # cross-fitted DML (Lasso + XGBoost)
├── 06_comparison.R                   # final comparison table + heatmap
├── 07_Robustness Check.R             # sparse-cell hypothesis: coarsening, bootstrap SE, seed sensitivity
├── 08_SparseCell Exclusion Check.R   # excludes the 10-row sparse cell, re-runs the comparison
├── 09_Simulation Check.R             # Monte Carlo check against a known true effect
├── Classic DiD vs DML.Rproj
├── data/
│   ├── public.dat                                    # original Card-Krueger survey data
│   ├── codebook                                      # variable definitions for public.dat
│   ├── read.me                                       # original data documentation
│   ├── survey1.nj / survey2.nj                       # original NJ survey instruments
│   ├── check.sas                                     # original SAS validation script
│   ├── fastfood_data.csv                             # cleaned full dataset
│   ├── estimation_sample.csv                         # primary complete-case sample
│   ├── estimation_sample_robustness_check.csv        # imputed/unified sample (missing = its own category)
│   ├── classic_did_summary.csv                       # primary DiD estimate (unified sample, HC1-robust)
│   ├── classic_did_summary_complete_case.csv         # reference DiD estimate (complete-case, classical SE)
│   └── dml_summary.csv                               # DML estimate summary
└── output/
    ├── figures/              # all plots (.png)
    └── tables/                # all regression/comparison tables (.csv), including table6-table10b from 07-09
```

## Reproducing the Results

Each script begins with `source(here::here("00_setup.R"))`, and the pipeline is meant to be run in numeric order:

```r
source("00_setup.R")
source("01_data_cleaning.R")
source("02_summary_tables.R")
source("03_classic_did_regression.R")
source("04_robustness check.R")
source("05_dml.R")
source("06_comparison.R")
source("07_Robustness Check.R")
source("08_SparseCell Exclusion Check.R")
source("09_Simulation Check.R")
```

Note that scripts 07-09 are computationally heavy: 07's bootstrap step alone re-fits the full DML cross-fitting procedure 500 times, and 09 does the same across 500 Monte Carlo replications. Each of these scripts has a `QUICK_TEST` flag near the top that can be set to `TRUE` to sanity-check the pipeline with far fewer iterations before committing to a full run.

Requirements: R with `readxl`, `dplyr`, `tidyr`, `lubridate`, `readr`, `stringr`, `here`, `ggplot2`, `broom`, `lmtest`, `sandwich`, `glmnet`, `xgboost`, `kableExtra` (auto-installed by `00_setup.R` if missing).

## References

The base classic-DiD replication (data construction and Tables 1-3) reproduces the original Card & Krueger methodology. The robustness check, DML estimation, final comparison framework, and the SE-reversal investigation (scripts 07-09) are our own extensions built on top of that replicated base.

- Card, D., and Krueger, A. B. (1994). Minimum Wages and Employment: A Case Study of the Fast-Food Industry in New Jersey and Pennsylvania. *American Economic Review*, 84(4), 772–793.
- Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W., and Robins, J. (2018). Double/Debiased Machine Learning for Treatment and Structural Parameters. *The Econometrics Journal*, 21(1), C1–C68.
- Robinson, P. M. (1988). Root-N-Consistent Semiparametric Regression. *Econometrica*, 56(4), 931–954.
