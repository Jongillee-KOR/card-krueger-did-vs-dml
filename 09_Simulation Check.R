#==================================================
# 9. Simulation Check: Monte Carlo with a Known True Effect
#==================================================
# 07-08 test the sparse-cell mechanism on the OBSERVED data -- a single real-
# world draw, so the "reversal" (and its disappearance under bootstrap SE)
# could in principle just be a fluke of this one sample. This script goes
# beyond observational results: it builds a synthetic DGP that reuses the
# REAL covariate/treatment joint distribution (so the sparse chain x
# ownership cell is preserved exactly, by resampling real rows rather than
# hand-tuning marginal probabilities) plus a real-data-calibrated nuisance
# signal, but imposes a KNOWN true treatment effect and fresh noise.
# Repeating this many times lets us directly measure each method's bias,
# variance, and RMSE against ground truth -- something the single real-data
# run in 03/05/06 cannot do (the true effect there is unknown). It also lets
# us check whether asymptotic HC1 SE is well-calibrated for each method
# (07/08 found it is NOT trustworthy on the real sample) by comparing the
# average reported HC1 SE against the actual Monte Carlo SD.
source(here::here("00_setup.R"))

QUICK_TEST <- FALSE
n_sims <- if (QUICK_TEST) 20 else 500
true_beta <- 12  # chosen, KNOWN ground truth -- deliberately inside the
# 11-14 range of the real point estimates (03/05/07), but
# fixed and known here so bias/variance can be measured
# against it. Not derived from the data.

#==================================================
# 9.1. Load real data & rebuild covariates (same as 03/05/07/08)
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
cat(sprintf("Loaded unified sample: N = %d\n", N))

#==================================================
# 9.2. Calibrate the DGP from the real data
#==================================================
# Fit the real unified-sample Classic DiD spec once, purely to extract a
# plausible nuisance signal (intercept + chain/ownership/region effects) and
# a residual SD for the simulation noise -- NOT used for inference about the
# real gap effect (that inference lives in 03/05/06). Here the gap
# coefficient from this fit is discarded and replaced by the chosen
# true_beta below.
calib_model <- lm(delta_emp ~ gap + chain_cat + co_owned_cat + centralj + northj + pa1, data = dml_df)
calib_coefs <- coef(calib_model)
calib_sigma <- summary(calib_model)$sigma

cat(sprintf("Calibration: N = %d, residual SD = %.4f (used as simulation noise SD)\n", N, calib_sigma))
cat(sprintf("Simulation ground truth: true_beta = %.2f (chosen, not estimated)\n", true_beta))

# Nuisance design matrix: same covariates as calib_model, minus gap. Column
# names line up with calib_coefs exactly since both come from the same data
# and default factor contrasts.
nuisance_design <- model.matrix(delta_emp ~ chain_cat + co_owned_cat + centralj + northj + pa1, data = dml_df)
nuisance_coefs  <- calib_coefs[colnames(nuisance_design)]
nuisance_signal <- as.numeric(nuisance_design %*% nuisance_coefs)  # per-row true intercept + covariate effect

#==================================================
# 9.3. Reusable DML cross-fit helper (same as 05/07/08)
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
# 9.4. Monte Carlo loop
#==================================================
cat(sprintf("\nRunning %d Monte Carlo replications...\n", n_sims))
X_vars <- model.matrix(~ chain_cat + co_owned_cat + centralj + northj + pa1 - 1, data = dml_df)

sim_results <- vector("list", n_sims)
for (s in 1:n_sims) {
  # Resample rows WITH replacement from the real sample. This preserves the
  # real (chain_cat, co_owned_cat, region, gap) joint distribution in every
  # replicate -- including the sparse Chain_4 x Co_Owned cell -- without
  # having to hand-tune marginal probabilities to reproduce it.
  boot_idx <- sample(seq_len(N), N, replace = TRUE)
  
  gap_sim      <- dml_df$gap[boot_idx]
  X_sim        <- X_vars[boot_idx, , drop = FALSE]
  nuisance_sim <- nuisance_signal[boot_idx]
  
  noise <- rnorm(N, mean = 0, sd = calib_sigma)
  Y_sim <- nuisance_sim + true_beta * gap_sim + noise
  
  sim_data <- dml_df[boot_idx, ] %>% mutate(delta_emp = Y_sim)
  
  # Classic DiD, HC1-robust
  m_classic <- lm(delta_emp ~ gap + chain_cat + co_owned_cat + centralj + northj + pa1, data = sim_data)
  m_classic_robust <- coeftest(m_classic, vcov = vcovHC(m_classic, type = "HC1"))
  classic_est <- m_classic_robust["gap", "Estimate"]
  classic_se  <- m_classic_robust["gap", "Std. Error"]
  
  # DML, same procedure as 05/07/08
  cf_sim  <- run_dml_crossfit(X_sim, Y_sim, gap_sim, seed = 2026 + s, K = 5)
  est_sim <- dml_estimate(cf_sim)
  
  sim_results[[s]] <- tibble(
    sim = s,
    classic_est = classic_est, classic_se = classic_se,
    lasso_est   = est_sim$lasso_est, lasso_se = est_sim$lasso_se,
    xgb_est     = est_sim$xgb_est,   xgb_se   = est_sim$xgb_se
  )
  
  if (s %% 50 == 0) cat(sprintf("  Simulation: %d / %d\n", s, n_sims))
}
sim_df <- bind_rows(sim_results)

#==================================================
# 9.5. Bias / variance / RMSE / coverage, per method
#==================================================
summarize_method <- function(est, se, true_val, label) {
  bias     <- mean(est) - true_val
  emp_sd   <- sd(est)                        # Monte Carlo "true" SE
  mean_se  <- mean(se)                       # average analytic (HC1) SE the method itself reports
  rmse     <- sqrt(mean((est - true_val)^2))
  ci_lower <- est - 1.96 * se
  ci_upper <- est + 1.96 * se
  coverage <- mean(true_val >= ci_lower & true_val <= ci_upper)
  tibble(method = label, bias = bias, empirical_sd = emp_sd,
         mean_reported_se = mean_se, rmse = rmse, coverage_95 = coverage)
}

table10_simulation <- bind_rows(
  summarize_method(sim_df$classic_est, sim_df$classic_se, true_beta, "Classic DiD"),
  summarize_method(sim_df$lasso_est,   sim_df$lasso_se,   true_beta, "DML-Lasso"),
  summarize_method(sim_df$xgb_est,     sim_df$xgb_se,     true_beta, "DML-XGBoost")
)

cat("\n=== Monte Carlo results (", n_sims, " replications, true_beta = ", true_beta, ") ===\n", sep = "")
print(table10_simulation)

#==================================================
# 9.6. Interpretation
#==================================================
classic_sd <- table10_simulation$empirical_sd[table10_simulation$method == "Classic DiD"]
lasso_sd   <- table10_simulation$empirical_sd[table10_simulation$method == "DML-Lasso"]
xgb_sd     <- table10_simulation$empirical_sd[table10_simulation$method == "DML-XGBoost"]

dml_lower_variance <- classic_sd > lasso_sd && classic_sd > xgb_sd

cat("\n--- Interpretation (simulation check) ---\n")
if (dml_lower_variance) {
  cat("Under a KNOWN true effect and the real sparse-cell covariate structure, DML's Monte Carlo SD is lower than Classic DiD's: contrary to 07/08's bootstrap-based finding, DML would be the more efficient estimator here in a fully controlled setting -- worth reconciling with the bootstrap results before drawing a final conclusion.\n")
} else {
  cat("Under simulation, Classic DiD's Monte Carlo SD is lower than (or equal to) DML's: this matches 07/08's bootstrap-based finding (the textbook prediction -- OLS more efficient in this low-dimensional, sparse-cell setting -- holds once SE is measured correctly), and confirms it wasn't just a feature of this one real sample.\n")
}

cat("\nBias check -- a method with |bias| noticeably larger than its empirical SD is trading meaningful bias for its variance reduction (relevant for judging whether any DML advantage is 'worth it'):\n")
print(table10_simulation %>% select(method, bias, empirical_sd))

cat("\nCalibration check -- mean_reported_se close to empirical_sd means the method's own HC1 SE is well-calibrated (near-nominal coverage); mean_reported_se noticeably below empirical_sd means the method understates its own uncertainty (this directly tests 07/08's finding that asymptotic HC1 SE is unreliable in this sample):\n")
print(table10_simulation %>% select(method, mean_reported_se, empirical_sd, coverage_95))

#==================================================
# 9.7. Plot & save
#==================================================
sim_long <- bind_rows(
  tibble(method = "Classic DiD", estimate = sim_df$classic_est),
  tibble(method = "DML-Lasso",   estimate = sim_df$lasso_est),
  tibble(method = "DML-XGBoost", estimate = sim_df$xgb_est)
)

p_sim <- ggplot(sim_long, aes(x = method, y = estimate, fill = method)) +
  geom_boxplot(alpha = 0.6, outlier.alpha = 0.3) +
  geom_hline(yintercept = true_beta, linetype = "dashed", color = "firebrick", linewidth = 1) +
  annotate("text", x = 0.6, y = true_beta, label = "True effect", color = "firebrick", vjust = -0.6, hjust = 0) +
  theme_classic() +
  labs(title = "Monte Carlo Distribution of Estimates vs Known True Effect",
       subtitle = sprintf("%d replications; real chain x ownership joint distribution preserved via resampling", n_sims),
       x = "Method", y = "Estimate") +
  theme(legend.position = "none")

print(p_sim)
ggsave(here("output", "figures", "simulation_check_boxplot.png"), p_sim, width = 8, height = 6, dpi = 300)

write_csv(table10_simulation, here("output", "tables", "table10_simulation_summary.csv"))
write_csv(sim_df,             here("output", "tables", "table10b_simulation_raw.csv"))

cat("\nSaved: output/tables/table10_simulation_summary.csv, table10b_simulation_raw.csv\n")
cat("Saved: output/figures/simulation_check_boxplot.png\n")