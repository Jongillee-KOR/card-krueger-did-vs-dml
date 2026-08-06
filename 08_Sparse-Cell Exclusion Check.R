#==================================================
# 8. Sparse-Cell Exclusion Check
#==================================================
# Original plan for this script split the sample into "sparse cell (n<20)"
# vs "sufficient sample (n>=20)" subgroups. Running that version showed the
# "sparse" subgroup was literally a single cell -- Chain_4 x Co_Owned, n=10
# (07's cell-count table: every other cell is >=20) -- so chain_cat and
# co_owned_cat are constant within it, and DML's cross-fitting couldn't be
# validated on ~5 obs/fold. See README Limitations for that design's outcome.
#
# This script replaces it with a lighter, more direct stress test: drop just
# those 10 rows from the FULL unified sample and rerun the Classic DiD vs DML
# comparison (asymptotic HC1 SE, per 03/05/06, AND bootstrap SE, per 07 --
# 07 showed the two can disagree) on the remaining N=360. Comparing against
# the full-sample (N=370) results from 07 shows directly how much that one
# cell was driving the SE pattern, without the subgroup-split's degeneracy
# problem.
source(here::here("00_setup.R"))

QUICK_TEST <- FALSE
n_boot_classic <- if (QUICK_TEST) 50 else 1000
n_boot_dml     <- if (QUICK_TEST) 20 else 500

#==================================================
# 8.1. Load data & rebuild covariates (same as 03/05/07)
#==================================================
data_path <- here("data", "estimation_sample_robustness_check.csv")
if (!file.exists(data_path)) {
  stop("Could not find ", data_path, ". Run 01_Data_cleaning.R then 04_Robustness_check.R first.")
}
est_df <- read_csv(data_path, show_col_types = FALSE)

dml_df <- est_df %>%
  filter(!is.na(delta_emp) & !is.na(gap)) %>%
  mutate(
    chain_cat = case_when(
      chain_was_na ~ "Missing",
      chain == 1 ~ "Chain_1",
      chain == 2 ~ "Chain_2",
      chain == 3 ~ "Chain_3",
      chain == 4 ~ "Chain_4",
      TRUE       ~ "Missing"
    ),
    co_owned_cat = case_when(
      co_owned_was_na ~ "Missing",
      co_owned == 1 ~ "Co_Owned",
      co_owned == 0 ~ "Franchise",
      TRUE          ~ "Missing"
    ),
    centralj = if_else(is.na(centralj), 0, centralj),
    northj   = if_else(is.na(northj), 0, northj),
    pa1      = if_else(is.na(pa1), 0, pa1)
  )

N_full <- nrow(dml_df)
cat(sprintf("Full unified sample: N = %d\n", N_full))

#==================================================
# 8.2. Identify & exclude the sparse cell
#==================================================
sparse_mask <- dml_df$chain_cat == "Chain_4" & dml_df$co_owned_cat == "Co_Owned"
cat(sprintf("Sparse cell (Chain_4 x Co_Owned) rows: %d\n", sum(sparse_mask)))

excl_df <- dml_df[!sparse_mask, ]
N_excl <- nrow(excl_df)
cat(sprintf("Excluded-sample N = %d (dropped %d rows from the full N = %d)\n",
            N_excl, N_full - N_excl, N_full))

# chain_cat still has all 5 levels (Chain_4 remains present via its
# Franchise rows, n=40) and co_owned_cat still has all 3 levels -- unlike the
# abandoned subgroup design, this does NOT create a degenerate design matrix.
cat("\nchain_cat levels remaining:", paste(sort(unique(excl_df$chain_cat)), collapse = ", "), "\n")
cat("co_owned_cat levels remaining:", paste(sort(unique(excl_df$co_owned_cat)), collapse = ", "), "\n")

#==================================================
# 8.3. Reusable DML cross-fitting helper (same as 05/07)
#==================================================
run_dml_crossfit <- function(X, Y, D, seed, K = 5) {
  set.seed(seed)
  n <- length(Y)
  folds <- sample(rep(1:K, length.out = n))
  
  D_tilde_lasso <- numeric(n); Y_tilde_lasso <- numeric(n)
  D_tilde_xgb   <- numeric(n); Y_tilde_xgb   <- numeric(n)
  
  for (k in 1:K) {
    idx_train <- which(folds != k)
    idx_test  <- which(folds == k)
    
    lasso_d <- cv.glmnet(X[idx_train, , drop = FALSE], D[idx_train], alpha = 1)
    D_tilde_lasso[idx_test] <- D[idx_test] -
      as.numeric(predict(lasso_d, newx = X[idx_test, , drop = FALSE], s = "lambda.min"))
    
    lasso_y <- cv.glmnet(X[idx_train, , drop = FALSE], Y[idx_train], alpha = 1)
    Y_tilde_lasso[idx_test] <- Y[idx_test] -
      as.numeric(predict(lasso_y, newx = X[idx_test, , drop = FALSE], s = "lambda.min"))
    
    xgb_d <- xgboost(data = X[idx_train, , drop = FALSE], label = D[idx_train],
                     objective = "reg:squarederror", nrounds = 100, max_depth = 3,
                     eta = 0.05, verbose = 0)
    D_tilde_xgb[idx_test] <- D[idx_test] - predict(xgb_d, X[idx_test, , drop = FALSE])
    
    xgb_y <- xgboost(data = X[idx_train, , drop = FALSE], label = Y[idx_train],
                     objective = "reg:squarederror", nrounds = 100, max_depth = 3,
                     eta = 0.05, verbose = 0)
    Y_tilde_xgb[idx_test] <- Y[idx_test] - predict(xgb_y, X[idx_test, , drop = FALSE])
  }
  
  list(D_tilde_lasso = D_tilde_lasso, Y_tilde_lasso = Y_tilde_lasso,
       D_tilde_xgb   = D_tilde_xgb,   Y_tilde_xgb   = Y_tilde_xgb)
}

dml_estimate <- function(cf) {
  m_lasso <- lm(cf$Y_tilde_lasso ~ cf$D_tilde_lasso - 1)
  m_xgb   <- lm(cf$Y_tilde_xgb   ~ cf$D_tilde_xgb   - 1)
  se_lasso <- coeftest(m_lasso, vcov = vcovHC(m_lasso, type = "HC1"))[1, "Std. Error"]
  se_xgb   <- coeftest(m_xgb,   vcov = vcovHC(m_xgb,   type = "HC1"))[1, "Std. Error"]
  list(lasso_est = unname(coef(m_lasso)[1]), lasso_se = unname(se_lasso),
       xgb_est   = unname(coef(m_xgb)[1]),   xgb_se   = unname(se_xgb))
}

#==================================================
# 8.4. Classic DiD & DML on the excluded sample, asymptotic HC1 SE
#==================================================
m_classic_excl <- lm(delta_emp ~ gap + chain_cat + co_owned_cat + centralj + northj + pa1, data = excl_df)
m_classic_excl_robust <- coeftest(m_classic_excl, vcov = vcovHC(m_classic_excl, type = "HC1"))
classic_excl_est <- m_classic_excl_robust["gap", "Estimate"]
classic_excl_se  <- m_classic_excl_robust["gap", "Std. Error"]

X_vars_excl <- model.matrix(~ chain_cat + co_owned_cat + centralj + northj + pa1 - 1, data = excl_df)
Y_excl <- excl_df$delta_emp
D_excl <- excl_df$gap

cf_excl  <- run_dml_crossfit(X_vars_excl, Y_excl, D_excl, seed = 2026, K = 5)
est_excl <- dml_estimate(cf_excl)

cat("\n=== Excluded sample (N =", N_excl, "), asymptotic HC1 SE ===\n", sep = "")
cat(sprintf("Classic DiD:  estimate = %.4f, SE = %.4f\n", classic_excl_est, classic_excl_se))
cat(sprintf("DML-Lasso:    estimate = %.4f, SE = %.4f\n", est_excl$lasso_est, est_excl$lasso_se))
cat(sprintf("DML-XGBoost:  estimate = %.4f, SE = %.4f\n", est_excl$xgb_est, est_excl$xgb_se))

#==================================================
# 8.5. Bootstrap SE on the excluded sample
#==================================================
cat(sprintf("\n--- Bootstrap SE on excluded sample (N=%d): Classic DiD %d reps | DML %d reps ---\n",
            N_excl, n_boot_classic, n_boot_dml))

boot_classic_est <- numeric(n_boot_classic)
for (b in 1:n_boot_classic) {
  boot_idx  <- sample(seq_len(N_excl), N_excl, replace = TRUE)
  boot_data <- excl_df[boot_idx, ]
  m_boot <- lm(delta_emp ~ gap + chain_cat + co_owned_cat + centralj + northj + pa1, data = boot_data)
  boot_classic_est[b] <- coef(m_boot)["gap"]
  if (b %% 200 == 0) cat(sprintf("  Classic DiD bootstrap: %d / %d\n", b, n_boot_classic))
}
bootstrap_classic_excl_se <- sd(boot_classic_est)

boot_lasso_est <- numeric(n_boot_dml)
boot_xgb_est   <- numeric(n_boot_dml)
for (b in 1:n_boot_dml) {
  boot_idx <- sample(seq_len(N_excl), N_excl, replace = TRUE)
  cf_boot  <- run_dml_crossfit(X_vars_excl[boot_idx, , drop = FALSE], Y_excl[boot_idx], D_excl[boot_idx],
                               seed = 4000 + b, K = 5)
  est_boot <- dml_estimate(cf_boot)
  boot_lasso_est[b] <- est_boot$lasso_est
  boot_xgb_est[b]   <- est_boot$xgb_est
  if (b %% 100 == 0) cat(sprintf("  DML bootstrap: %d / %d\n", b, n_boot_dml))
}
bootstrap_lasso_excl_se <- sd(boot_lasso_est)
bootstrap_xgb_excl_se   <- sd(boot_xgb_est)

cat("\n=== Excluded sample, bootstrap SE ===\n")
cat(sprintf("Classic DiD:  bootstrap SE = %.4f\n", bootstrap_classic_excl_se))
cat(sprintf("DML-Lasso:    bootstrap SE = %.4f\n", bootstrap_lasso_excl_se))
cat(sprintf("DML-XGBoost:  bootstrap SE = %.4f\n", bootstrap_xgb_excl_se))

#==================================================
# 8.6. Load 07's full-sample (N=370) results for comparison
#==================================================
classic_path <- here("data", "classic_did_summary.csv")
dml_path     <- here("data", "dml_summary.csv")
bootstrap_path <- here("output", "tables", "table7_bootstrap_se.csv")
if (!file.exists(classic_path))   stop("Could not find ", classic_path, ". Run 03_Classic_Did_Regression.R first.")
if (!file.exists(dml_path))       stop("Could not find ", dml_path,     ". Run 05_DML.R first.")
if (!file.exists(bootstrap_path)) stop("Could not find ", bootstrap_path, ". Run 07_Robustness_checks.R first.")

baseline_classic  <- read_csv(classic_path, show_col_types = FALSE)
baseline_dml      <- read_csv(dml_path,     show_col_types = FALSE)
baseline_bootstrap <- read_csv(bootstrap_path, show_col_types = FALSE)

full_classic_hc1  <- baseline_classic$std_error[1]
full_lasso_hc1    <- baseline_dml$std_error[baseline_dml$method == "DML (Lasso)"]
full_xgb_hc1      <- baseline_dml$std_error[baseline_dml$method == "DML (XGBoost)"]

full_classic_boot <- baseline_bootstrap$bootstrap_se[baseline_bootstrap$method == "Classic DiD"]
full_lasso_boot   <- baseline_bootstrap$bootstrap_se[baseline_bootstrap$method == "DML-Lasso"]
full_xgb_boot     <- baseline_bootstrap$bootstrap_se[baseline_bootstrap$method == "DML-XGBoost"]

#==================================================
# 8.7. Combined comparison table
#==================================================
table9_exclusion <- tibble(
  sample  = c(rep("Full (N=370)", 3), rep("Sparse cell excluded (N=360)", 3)),
  method  = rep(c("Classic DiD", "DML-Lasso", "DML-XGBoost"), 2),
  asymptotic_hc1_se = c(full_classic_hc1, full_lasso_hc1, full_xgb_hc1,
                        classic_excl_se, est_excl$lasso_se, est_excl$xgb_se),
  bootstrap_se       = c(full_classic_boot, full_lasso_boot, full_xgb_boot,
                         bootstrap_classic_excl_se, bootstrap_lasso_excl_se, bootstrap_xgb_excl_se)
)

cat("\n=== Full sample vs sparse-cell-excluded sample: HC1 SE and bootstrap SE ===\n")
print(table9_exclusion)

#==================================================
# 8.8. Interpretation
#==================================================
cat("\n--- Interpretation (sparse-cell exclusion) ---\n")

hc1_reversal_full <- full_classic_hc1 > full_lasso_hc1 && full_classic_hc1 > full_xgb_hc1
hc1_reversal_excl <- classic_excl_se > est_excl$lasso_se && classic_excl_se > est_excl$xgb_se
boot_reversal_full <- full_classic_boot > full_lasso_boot && full_classic_boot > full_xgb_boot
boot_reversal_excl <- bootstrap_classic_excl_se > bootstrap_lasso_excl_se && bootstrap_classic_excl_se > bootstrap_xgb_excl_se

cat(sprintf("Asymptotic HC1 reversal (Classic DiD SE largest): full sample = %s, excluded sample = %s\n",
            hc1_reversal_full, hc1_reversal_excl))
cat(sprintf("Bootstrap reversal (Classic DiD SE largest):      full sample = %s, excluded sample = %s\n",
            boot_reversal_full, boot_reversal_excl))

if (hc1_reversal_full != hc1_reversal_excl || boot_reversal_full != boot_reversal_excl) {
  cat("\nDropping just the 10-row sparse cell CHANGES the reversal pattern (on at least one SE estimator): this single cell has a disproportionate influence relative to its 10/370 (~2.7%) share of the sample -- meaningful support for the sparse-cell explanation.\n")
} else {
  cat("\nDropping the 10-row sparse cell does NOT change the reversal pattern on either SE estimator: the pattern (whichever direction it runs) is not driven by this one cell specifically -- it reflects something about the sample/estimators more broadly, consistent with 07's finding that asymptotic HC1 SE itself is the more likely culprit rather than any single sparse cell.\n")
}

#==================================================
# 8.9. Plot & save
#==================================================
table9_exclusion_long <- table9_exclusion %>%
  pivot_longer(cols = c(asymptotic_hc1_se, bootstrap_se), names_to = "se_type", values_to = "se") %>%
  mutate(se_type = if_else(se_type == "asymptotic_hc1_se", "Asymptotic HC1", "Bootstrap"))

p_exclusion <- ggplot(table9_exclusion_long, aes(x = method, y = se, fill = sample)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  facet_wrap(~ se_type) +
  theme_classic() +
  labs(title = "SE by Method: Full Sample vs Sparse-Cell-Excluded Sample",
       subtitle = "Sparse cell = Chain_4 x Co_Owned, n = 10 (dropped in the right-hand bars)",
       x = "Method", y = "Standard Error", fill = "Sample") +
  theme(legend.position = "bottom")

print(p_exclusion)
ggsave(here("output", "figures", "sparse_cell_exclusion_comparison.png"), p_exclusion, width = 9, height = 6, dpi = 300)
write_csv(table9_exclusion, here("output", "tables", "table9_sparse_cell_exclusion.csv"))

cat("\nSaved: output/tables/table9_sparse_cell_exclusion.csv\n")
cat("Saved: output/figures/sparse_cell_exclusion_comparison.png\n")