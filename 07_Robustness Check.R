#==================================================
# 7. Robustness Checks: Testing the Sparse-Cell Hypothesis
#==================================================
# Working hypothesis (from 03/05/06): the SE reversal on the unified sample
# (Classic DiD SE 7.13 > DML-Lasso 6.98 > DML-XGBoost 6.85) is driven by the
# 8-cell chain x ownership interaction implied by OLS's dummy encoding (min
# cell n=10, Chain_4 x Co_Owned), which DML's regularized first stage
# handles more gracefully than OLS. This script runs three checks:
#   7.4 Coarser categorical: collapse chain to 2 levels, does the reversal disappear?
#   7.5 Bootstrap SE: does the reversal survive a resampling-based SE instead of asymptotic HC1?
#   7.6 Seed sensitivity: is the DML side of the reversal robust to the 5-fold split draw?
source(here::here("00_setup.R"))

# QUICK_TEST: set TRUE to sanity-check the pipeline with far fewer iterations
# before committing to the full 500-1000 iteration run (the DML bootstrap in
# 7.5 re-runs the entire cross-fitting procedure per replicate and is slow).
QUICK_TEST <- FALSE

n_boot_classic <- if (QUICK_TEST) 50 else 1000
n_boot_dml     <- if (QUICK_TEST) 20 else 500
n_seeds        <- if (QUICK_TEST) 5  else 30

#==================================================
# 7.1. Load data & rebuild covariates (same as 03 section 3.7 / 05 sections 5.1-5.2)
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

N <- nrow(dml_df)
Y <- dml_df$delta_emp
D <- dml_df$gap
X_vars <- model.matrix(~ chain_cat + co_owned_cat + centralj + northj + pa1 - 1, data = dml_df)

cat(sprintf("Loaded unified sample: N = %d\n", N))

#==================================================
# 7.2. Load baseline (5-level chain) results for reference
#==================================================
classic_path <- here("data", "classic_did_summary.csv")
dml_path     <- here("data", "dml_summary.csv")
if (!file.exists(classic_path)) stop("Could not find ", classic_path, ". Run 03_Classic_Did_Regression.R first.")
if (!file.exists(dml_path))     stop("Could not find ", dml_path,     ". Run 05_DML.R first.")

baseline_classic <- read_csv(classic_path, show_col_types = FALSE)
baseline_dml     <- read_csv(dml_path,     show_col_types = FALSE)

baseline_classic_se <- baseline_classic$std_error[1]
baseline_lasso_se   <- baseline_dml$std_error[baseline_dml$method == "DML (Lasso)"]
baseline_xgb_se     <- baseline_dml$std_error[baseline_dml$method == "DML (XGBoost)"]

cat("\n=== Baseline (5-level chain, N =", N, ") ===\n", sep = "")
cat(sprintf("Classic DiD SE: %.4f | DML-Lasso SE: %.4f | DML-XGBoost SE: %.4f\n",
            baseline_classic_se, baseline_lasso_se, baseline_xgb_se))
cat("(Baseline shows the reversal: Classic DiD SE is the largest of the three.)\n")

#==================================================
# 7.3. Reusable DML cross-fitting helper (used in 7.4, 7.5, 7.6)
#==================================================
# Same 5-fold cross-fitted Lasso & XGBoost nuisance procedure as 05_DML.R,
# parameterized by covariate matrix / outcome / treatment / seed, so it can be
# re-run under different encodings, resampled data, and fold-assignment seeds
# without duplicating the fitting logic.
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

# Given cross-fitted residuals, returns point estimates and HC1-robust SEs
# for both learners via residual-on-residual regression (same as 05 section 5.4).
dml_estimate <- function(cf) {
  m_lasso <- lm(cf$Y_tilde_lasso ~ cf$D_tilde_lasso - 1)
  m_xgb   <- lm(cf$Y_tilde_xgb   ~ cf$D_tilde_xgb   - 1)
  se_lasso <- coeftest(m_lasso, vcov = vcovHC(m_lasso, type = "HC1"))[1, "Std. Error"]
  se_xgb   <- coeftest(m_xgb,   vcov = vcovHC(m_xgb,   type = "HC1"))[1, "Std. Error"]
  list(lasso_est = unname(coef(m_lasso)[1]), lasso_se = unname(se_lasso),
       xgb_est   = unname(coef(m_xgb)[1]),   xgb_se   = unname(se_xgb))
}

#==================================================
# 7.4. Coarser categorical check: 5-level chain -> 2-level (major chain vs other)
#==================================================
cat("\n--- 7.4 Coarser categorical check ---\n")

# 7.4.1 Cell counts before coarsening (confirms the sparse-cell structure)
cat("\nCell counts, chain_cat x co_owned_cat (5-level, before coarsening):\n")
print(table(dml_df$chain_cat, dml_df$co_owned_cat))

# 7.4.2 Collapse chain into 2 levels: the most frequent actual chain
# ("Major_Chain") vs everything else (minor chains + Missing). This folds the
# sparse Chain_4 x Co_Owned cell into a much larger "Other" cell.
chain_counts <- table(dml_df$chain_cat)
chain_counts_known <- chain_counts[names(chain_counts) != "Missing"]
major_chain_name <- names(which.max(chain_counts_known))
cat(sprintf("\nMost frequent chain (used as 'Major_Chain'): %s (n = %d)\n",
            major_chain_name, chain_counts_known[major_chain_name]))

dml_df <- dml_df %>%
  mutate(chain_cat2 = if_else(chain_cat == major_chain_name, "Major_Chain", "Other"))

cat("\nCell counts, chain_cat2 x co_owned_cat (2-level, after coarsening):\n")
print(table(dml_df$chain_cat2, dml_df$co_owned_cat))

# 7.4.3 Classic DiD, coarser encoding, HC1-robust SE
model7_classic_coarser <- lm(delta_emp ~ gap + chain_cat2 + co_owned_cat + centralj + northj + pa1,
                             data = dml_df)
model7_classic_coarser_robust <- coeftest(model7_classic_coarser,
                                          vcov = vcovHC(model7_classic_coarser, type = "HC1"))
coarser_classic_est <- model7_classic_coarser_robust["gap", "Estimate"]
coarser_classic_se  <- model7_classic_coarser_robust["gap", "Std. Error"]

# 7.4.4 DML, coarser encoding, same cross-fitting procedure as 05, seed = 2026
X_vars2 <- model.matrix(~ chain_cat2 + co_owned_cat + centralj + northj + pa1 - 1, data = dml_df)

cf_coarser  <- run_dml_crossfit(X_vars2, Y, D, seed = 2026, K = 5)
est_coarser <- dml_estimate(cf_coarser)

cat("\n=== Coarser categorical (2-level chain, N =", N, ") ===\n", sep = "")
cat(sprintf("Classic DiD:  estimate = %.4f, SE = %.4f\n", coarser_classic_est, coarser_classic_se))
cat(sprintf("DML-Lasso:    estimate = %.4f, SE = %.4f\n", est_coarser$lasso_est, est_coarser$lasso_se))
cat(sprintf("DML-XGBoost:  estimate = %.4f, SE = %.4f\n", est_coarser$xgb_est, est_coarser$xgb_se))

table6_coarser <- tibble(
  encoding = c(rep("5-level (baseline)", 3), rep("2-level (coarser)", 3)),
  method   = rep(c("Classic DiD", "DML-Lasso", "DML-XGBoost"), 2),
  estimate = c(baseline_classic$estimate[1],
               baseline_dml$estimate[baseline_dml$method == "DML (Lasso)"],
               baseline_dml$estimate[baseline_dml$method == "DML (XGBoost)"],
               coarser_classic_est, est_coarser$lasso_est, est_coarser$xgb_est),
  std_error = c(baseline_classic_se, baseline_lasso_se, baseline_xgb_se,
                coarser_classic_se, est_coarser$lasso_se, est_coarser$xgb_se)
)

cat("\n--- Interpretation (coarser categorical check) ---\n")
if (coarser_classic_se < est_coarser$lasso_se && coarser_classic_se < est_coarser$xgb_se) {
  cat("Coarsening restores the original ordering (Classic DiD SE < DML SE): supports the sparse-cell explanation.\n")
} else if (coarser_classic_se < baseline_classic_se) {
  cat("Classic DiD SE fell after coarsening but did not fall below the DML SEs: partial support only.\n")
} else {
  cat("Classic DiD SE did not fall after coarsening: this check does NOT support the sparse-cell explanation.\n")
}

p_coarser <- ggplot(table6_coarser, aes(x = method, y = estimate, color = encoding)) +
  geom_point(position = position_dodge(width = 0.3), size = 4) +
  geom_errorbar(aes(ymin = estimate - 1.96 * std_error, ymax = estimate + 1.96 * std_error),
                position = position_dodge(width = 0.3), width = 0.15, linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  theme_classic() +
  labs(title = "5-Level vs 2-Level Chain Encoding: Classic DiD vs DML",
       subtitle = "If the sparse-cell explanation holds, Classic DiD's error bars should shrink under 2-level encoding",
       x = "Method", y = "Estimate & 95% CI", color = "Chain Encoding")

print(p_coarser)
ggsave(here("output", "figures", "coarser_categorical_comparison.png"), p_coarser, width = 8, height = 6, dpi = 300)
write_csv(table6_coarser, here("output", "tables", "table6_coarser_categorical.csv"))

#==================================================
# 7.5. Bootstrap SE: does the reversal survive a resampling-based SE?
#==================================================
cat("\n--- 7.5 Bootstrap SE ---\n")
cat(sprintf("Classic DiD: %d bootstrap resamples | DML: %d bootstrap resamples (5-level baseline encoding)\n",
            n_boot_classic, n_boot_dml))

# 7.5.1 Classic DiD bootstrap (original 5-level encoding, unified sample)
boot_classic_est <- numeric(n_boot_classic)
for (b in 1:n_boot_classic) {
  boot_idx  <- sample(seq_len(N), N, replace = TRUE)
  boot_data <- dml_df[boot_idx, ]
  m_boot <- lm(delta_emp ~ gap + chain_cat + co_owned_cat + centralj + northj + pa1, data = boot_data)
  boot_classic_est[b] <- coef(m_boot)["gap"]
  if (b %% 100 == 0) cat(sprintf("  Classic DiD bootstrap: %d / %d\n", b, n_boot_classic))
}
bootstrap_classic_se <- sd(boot_classic_est)

# 7.5.2 DML bootstrap (original 5-level encoding). The full cross-fitting
# pipeline is re-run on each resample, since cross-fitting (fold assignment +
# nuisance model refit) is part of the estimator, not just the final step.
boot_lasso_est <- numeric(n_boot_dml)
boot_xgb_est   <- numeric(n_boot_dml)
for (b in 1:n_boot_dml) {
  boot_idx <- sample(seq_len(N), N, replace = TRUE)
  cf_boot  <- run_dml_crossfit(X_vars[boot_idx, , drop = FALSE], Y[boot_idx], D[boot_idx],
                               seed = 1000 + b, K = 5)
  est_boot <- dml_estimate(cf_boot)
  boot_lasso_est[b] <- est_boot$lasso_est
  boot_xgb_est[b]   <- est_boot$xgb_est
  if (b %% 50 == 0) cat(sprintf("  DML bootstrap: %d / %d\n", b, n_boot_dml))
}
bootstrap_lasso_se <- sd(boot_lasso_est)
bootstrap_xgb_se   <- sd(boot_xgb_est)

table7_bootstrap <- tibble(
  method            = c("Classic DiD", "DML-Lasso", "DML-XGBoost"),
  asymptotic_hc1_se = c(baseline_classic_se, baseline_lasso_se, baseline_xgb_se),
  bootstrap_se      = c(bootstrap_classic_se, bootstrap_lasso_se, bootstrap_xgb_se),
  n_boot            = c(n_boot_classic, n_boot_dml, n_boot_dml)
)

cat("\n=== Asymptotic HC1 SE vs Bootstrap SE ===\n")
print(table7_bootstrap)

cat("\n--- Interpretation (bootstrap SE) ---\n")
if (bootstrap_classic_se > bootstrap_lasso_se && bootstrap_classic_se > bootstrap_xgb_se) {
  cat("The reversal reproduces under bootstrap resampling: Classic DiD's bootstrap SE is still the largest.\n")
} else {
  cat("The reversal does NOT reproduce under bootstrap resampling -- the asymptotic HC1 estimate may be unreliable in this small, sparse-cell sample.\n")
}

boot_long <- bind_rows(
  tibble(method = "Classic DiD", estimate = boot_classic_est),
  tibble(method = "DML-Lasso",   estimate = boot_lasso_est),
  tibble(method = "DML-XGBoost", estimate = boot_xgb_est)
)

p_bootstrap <- ggplot(boot_long, aes(x = method, y = estimate, fill = method)) +
  geom_boxplot(alpha = 0.6, outlier.alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  theme_classic() +
  labs(title = "Bootstrap Distribution of the Gap Coefficient",
       subtitle = sprintf("Classic DiD: %d resamples | DML: %d resamples", n_boot_classic, n_boot_dml),
       x = "Method", y = "Bootstrap Estimate") +
  theme(legend.position = "none")

print(p_bootstrap)
ggsave(here("output", "figures", "bootstrap_se_comparison.png"), p_bootstrap, width = 8, height = 6, dpi = 300)
write_csv(table7_bootstrap, here("output", "tables", "table7_bootstrap_se.csv"))

#==================================================
# 7.6. Seed sensitivity: does the reversal depend on one lucky fold split?
#==================================================
cat("\n--- 7.6 Seed sensitivity ---\n")
seeds <- seq_len(n_seeds) * 7 + 100  # arbitrary spread of seeds, distinct from seed = 2026 used in 05

seed_results <- vector("list", n_seeds)
for (i in seq_along(seeds)) {
  cf_seed  <- run_dml_crossfit(X_vars, Y, D, seed = seeds[i], K = 5)
  est_seed <- dml_estimate(cf_seed)
  seed_results[[i]] <- tibble(seed = seeds[i],
                              lasso_est = est_seed$lasso_est, lasso_se = est_seed$lasso_se,
                              xgb_est   = est_seed$xgb_est,   xgb_se   = est_seed$xgb_se)
  if (i %% 5 == 0) cat(sprintf("  Seed sensitivity: %d / %d\n", i, n_seeds))
}
seed_df <- bind_rows(seed_results)

cat("\n=== DML estimate/SE across", n_seeds, "fold-assignment seeds ===\n")
cat(sprintf("Lasso SE:    mean = %.4f, sd = %.4f, range = [%.4f, %.4f]\n",
            mean(seed_df$lasso_se), sd(seed_df$lasso_se), min(seed_df$lasso_se), max(seed_df$lasso_se)))
cat(sprintf("XGBoost SE:  mean = %.4f, sd = %.4f, range = [%.4f, %.4f]\n",
            mean(seed_df$xgb_se), sd(seed_df$xgb_se), min(seed_df$xgb_se), max(seed_df$xgb_se)))
cat(sprintf("Classic DiD SE (fixed, no fold randomness): %.4f\n", baseline_classic_se))

n_lasso_exceeds <- sum(seed_df$lasso_se > baseline_classic_se)
n_xgb_exceeds   <- sum(seed_df$xgb_se   > baseline_classic_se)

cat("\n--- Interpretation (seed sensitivity) ---\n")
if (n_lasso_exceeds == 0 && n_xgb_exceeds == 0) {
  cat("DML's SE stays below Classic DiD's SE across all seeds tested: the reversal is not an artifact of one lucky fold split.\n")
} else {
  cat(sprintf("DML's SE exceeded Classic DiD's SE in %d/%d seeds (Lasso) and %d/%d seeds (XGBoost): the reversal is partly seed-dependent.\n",
              n_lasso_exceeds, n_seeds, n_xgb_exceeds, n_seeds))
}

write_csv(seed_df, here("output", "tables", "table8_seed_sensitivity.csv"))

seed_long <- seed_df %>%
  select(seed, lasso_se, xgb_se) %>%
  pivot_longer(cols = c(lasso_se, xgb_se), names_to = "method", values_to = "se") %>%
  mutate(method = if_else(method == "lasso_se", "DML-Lasso", "DML-XGBoost"))

p_seed <- ggplot(seed_long, aes(x = method, y = se, fill = method)) +
  geom_boxplot(alpha = 0.6) +
  geom_hline(yintercept = baseline_classic_se, linetype = "dashed", color = "firebrick", linewidth = 1) +
  annotate("text", x = 1.5, y = baseline_classic_se, label = "Classic DiD SE",
           vjust = -0.5, color = "firebrick") +
  theme_classic() +
  labs(title = "DML SE Across Fold-Assignment Seeds",
       subtitle = sprintf("%d seeds; dashed line = Classic DiD's (fixed) SE", n_seeds),
       x = "Method", y = "Standard Error") +
  theme(legend.position = "none")

print(p_seed)
ggsave(here("output", "figures", "seed_sensitivity_boxplot.png"), p_seed, width = 8, height = 6, dpi = 300)

#==================================================
# 7.7. Combined summary
#==================================================
cat("\n=================================================\n")
cat("07_Robustness_checks.R summary\n")
cat("=================================================\n")
cat(sprintf("(1) Coarser categorical: Classic DiD SE %.4f -> %.4f (5-level baseline -> 2-level)\n",
            baseline_classic_se, coarser_classic_se))
cat(sprintf("(2) Bootstrap SE: Classic DiD %.4f, Lasso %.4f, XGBoost %.4f (vs asymptotic HC1 %.4f / %.4f / %.4f)\n",
            bootstrap_classic_se, bootstrap_lasso_se, bootstrap_xgb_se,
            baseline_classic_se, baseline_lasso_se, baseline_xgb_se))
cat(sprintf("(3) Seed sensitivity: Lasso SE range [%.4f, %.4f], XGBoost SE range [%.4f, %.4f], Classic DiD SE (fixed) %.4f\n",
            min(seed_df$lasso_se), max(seed_df$lasso_se), min(seed_df$xgb_se), max(seed_df$xgb_se), baseline_classic_se))
cat("\nSaved: output/tables/table6_coarser_categorical.csv, table7_bootstrap_se.csv, table8_seed_sensitivity.csv\n")
cat("Saved: output/figures/coarser_categorical_comparison.png, bootstrap_se_comparison.png, seed_sensitivity_boxplot.png\n")