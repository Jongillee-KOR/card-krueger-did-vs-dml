#==================================================
# 5. Double Machine Learning (Cross-Fitted)
#==================================================
source(here::here("00_setup.R"))

#==================================================
# 5.1. Load data
#==================================================
# Note: Uses the imputed/robustness sample (04), not the complete-case sample (03).
# Rationale: Retains explicit missingness flags (e.g., _was_na) so ML models 
# encode missingness as categories rather than dropping rows like lm() does, 
# allowing us to utilize more observations than standard regressions.
data_path <- here("data", "estimation_sample_robustness_check.csv")
if (!file.exists(data_path)) {
  stop("Could not find ", data_path, ". Run 04_sensitivity_imputation.R first.")
}
est_df <- read_csv(data_path, show_col_types = FALSE)

#==================================================
# 5.2. Build covariate matrix, outcome, and treatment
#==================================================
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
    centralj = if_else(is.na(centralj), 0, centralj), # Missing region dummies treated as 0
    northj   = if_else(is.na(northj), 0, northj),
    pa1      = if_else(is.na(pa1), 0, pa1)
  )

# X (covariates), Y (outcome), D (treatment: initial wage gap)
X_vars <- model.matrix(~ chain_cat + co_owned_cat + centralj + northj + pa1 - 1, data = dml_df)
Y <- dml_df$delta_emp
D <- dml_df$gap
N <- length(Y)

cat(sprintf("DML sample size: %d\n", N))
cat("Feature matrix columns:", paste(colnames(X_vars), collapse = ", "), "\n")

#==================================================
# 5.3. Cross-fitted nuisance function estimation (Lasso & XGBoost)
#==================================================
set.seed(2026)
K <- 5
folds <- sample(rep(1:K, length.out = N))

D_tilde_lasso <- numeric(N)
Y_tilde_lasso <- numeric(N)
D_tilde_xgb   <- numeric(N)
Y_tilde_xgb   <- numeric(N)

for (k in 1:K) {
  idx_train <- which(folds != k)
  idx_test  <- which(folds == k)
  
  # Lasso nuisance models
  lasso_d <- cv.glmnet(X_vars[idx_train, ], D[idx_train], alpha = 1)
  D_tilde_lasso[idx_test] <- D[idx_test] - as.numeric(predict(lasso_d, newx = X_vars[idx_test, ], s = "lambda.min"))
  
  lasso_y <- cv.glmnet(X_vars[idx_train, ], Y[idx_train], alpha = 1)
  Y_tilde_lasso[idx_test] <- Y[idx_test] - as.numeric(predict(lasso_y, newx = X_vars[idx_test, ], s = "lambda.min"))
  
  # XGBoost nuisance models
  xgb_d <- xgboost(data = X_vars[idx_train, ], label = D[idx_train],
                   objective = "reg:squarederror", nrounds = 100, max_depth = 3, eta = 0.05, verbose = 0)
  D_tilde_xgb[idx_test] <- D[idx_test] - predict(xgb_d, X_vars[idx_test, ])
  
  xgb_y <- xgboost(data = X_vars[idx_train, ], label = Y[idx_train],
                   objective = "reg:squarederror", nrounds = 100, max_depth = 3, eta = 0.05, verbose = 0)
  Y_tilde_xgb[idx_test] <- Y[idx_test] - predict(xgb_y, X_vars[idx_test, ])
}

cat("Cross-fitted nuisance function estimation completed.\n")

#==================================================
# 5.4. Final estimation: residual-on-residual regression
#==================================================
dml_lasso_model <- lm(Y_tilde_lasso ~ D_tilde_lasso - 1)
dml_xgb_model   <- lm(Y_tilde_xgb ~ D_tilde_xgb - 1)

res_lasso <- coeftest(dml_lasso_model, vcov = vcovHC(dml_lasso_model, type = "HC1"))
res_xgb   <- coeftest(dml_xgb_model, vcov = vcovHC(dml_xgb_model, type = "HC1"))

cat("\n=== DML (Lasso, cross-fitted) ===\n")
print(res_lasso)
cat("\n=== DML (XGBoost, cross-fitted) ===\n")
print(res_xgb)

#==================================================
# 5.5. Residual scatter plots
#==================================================
res_df <- tibble(
  D_tilde_lasso = D_tilde_lasso, Y_tilde_lasso = Y_tilde_lasso,
  D_tilde_xgb   = D_tilde_xgb,   Y_tilde_xgb   = Y_tilde_xgb
)

fig_res_lasso <- ggplot(res_df, aes(x = D_tilde_lasso, y = Y_tilde_lasso)) +
  geom_point(alpha = 0.3, color = "steelblue") +
  geom_smooth(method = "lm", formula = y ~ x, color = "firebrick", se = FALSE, linewidth = 1) +
  theme_classic() +
  labs(title = "DML-Lasso (Cross-Fitted)",
       subtitle = paste("Orthogonalized Treatment vs Outcome (Slope =", round(coef(dml_lasso_model), 4), ")"),
       x = "Treatment Residual (D_tilde)", y = "Outcome Residual (Y_tilde)")

fig_res_xgb <- ggplot(res_df, aes(x = D_tilde_xgb, y = Y_tilde_xgb)) +
  geom_point(alpha = 0.3, color = "darkgreen") +
  geom_smooth(method = "lm", formula = y ~ x, color = "firebrick", se = FALSE, linewidth = 1) +
  theme_classic() +
  labs(title = "DML-XGBoost (Cross-Fitted)",
       subtitle = paste("Orthogonalized Treatment vs Outcome (Slope =", round(coef(dml_xgb_model), 4), ")"),
       x = "Treatment Residual (D_tilde)", y = "Outcome Residual (Y_tilde)")

print(fig_res_lasso)
print(fig_res_xgb)
ggsave(here("output", "figures", "dml_lasso_residuals.png"), fig_res_lasso, width = 6, height = 4.5, dpi = 300)
ggsave(here("output", "figures", "dml_xgboost_residuals.png"), fig_res_xgb, width = 6, height = 4.5, dpi = 300)

#==================================================
# 5.6. Save outputs
#==================================================
dml_summary <- tibble(
  method    = c("DML (Lasso)", "DML (XGBoost)"),
  estimate  = c(res_lasso["D_tilde_lasso", "Estimate"], res_xgb["D_tilde_xgb", "Estimate"]),
  std_error = c(res_lasso["D_tilde_lasso", "Std. Error"], res_xgb["D_tilde_xgb", "Std. Error"])
)

write_csv(dml_summary, here("data", "dml_summary.csv"))
cat("\nSaved: data/dml_summary.csv, output/figures/dml_lasso_residuals.png, output/figures/dml_xgboost_residuals.png\n")