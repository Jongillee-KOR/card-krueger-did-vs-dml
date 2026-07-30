#==================================================
# 4. Robustness check: Missing Values Imputed as 0
#==================================================
source(here::here("00_setup.R"))

#=============================================================
# 4.1. Examine missing values in the original cleaned data
#=============================================================
df <- read_csv(here("data", "fastfood_data.csv"), show_col_types = FALSE)

print.data.frame(
  df %>% filter(is.na(empft) == TRUE & is.na(nmgrs) == FALSE)
)

#===========================================================
# 4.2. Alternative approach: impute all missing values as 0
#===========================================================
df_imputed <- df %>%
  mutate(
    co_owned_was_na = is.na(co_owned),
    chain_was_na    = is.na(chain)
  ) %>%
  replace(is.na(.), 0)

df_imputed <- df_imputed %>%
  mutate(
    fte  = empft + nmgrs + 0.5 * emppt,     # recalculate baseline FTE with imputed values
    fte2 = empft2 + nmgrs2 + 0.5 * emppt2   # recalculate follow-up FTE with imputed values
  ) %>%
  filter(fte > 0, fte2 > 0, wage_st > 0, wage_st2 > 0)

est_df_imputed <- df_imputed %>%
  ungroup() %>%
  mutate(
    delta_emp = fte2 - fte,
    gap = if_else(state == 1 & wage_st <= 5.05, (5.05 - wage_st) / wage_st, 0),
    chain1 = if_else(chain == 1, 1, 0),
    chain2 = if_else(chain == 2, 1, 0),
    chain3 = if_else(chain == 3, 1, 0),
    chain4 = if_else(chain == 4, 1, 0)
  )

cat(sprintf("Imputed-sample rows: %d (primary estimation_sample.csv had a different row count -- see 01_data_cleaning.R)\n",
            nrow(est_df_imputed)))

#==================================================
# 4.3. Re-run the regressions on the imputed sample
#==================================================
model1_imp <- lm(delta_emp ~ state, data = est_df_imputed)
model2_imp <- lm(delta_emp ~ state + co_owned + chain2 + chain3 + chain4, data = est_df_imputed)
model3_imp <- lm(delta_emp ~ gap, data = est_df_imputed)
model4_imp <- lm(delta_emp ~ gap + co_owned + chain2 + chain3 + chain4, data = est_df_imputed)

model1_imp_coeffs <- summary(model1_imp)$coefficients
model2_imp_coeffs <- summary(model2_imp)$coefficients
model3_imp_coeffs <- summary(model3_imp)$coefficients
model4_imp_coeffs <- summary(model4_imp)$coefficients

model1_imp_sigma <- summary(model1_imp)$sigma
model2_imp_sigma <- summary(model2_imp)$sigma
model3_imp_sigma <- summary(model3_imp)$sigma
model4_imp_sigma <- summary(model4_imp)$sigma

table4_imputed <- data.frame(
  Independent_Var = c("New Jersey Dummy", "Initial Wage Gap", "Controls for Chain and Ownership",
                      "Standard Error of Regression", "Number of Stores in Sample"),
  Model1_Coeff = c(round(model1_imp_coeffs["state", "Estimate"], 2), "-", "No",
                   round(model1_imp_sigma, 2), nobs(model1_imp)),
  Model2_Coeff = c(round(model2_imp_coeffs["state", "Estimate"], 2), "-", "Yes",
                   round(model2_imp_sigma, 2), nobs(model2_imp)),
  Model3_Coeff = c("-", round(model3_imp_coeffs["gap", "Estimate"], 2), "No",
                   round(model3_imp_sigma, 2), nobs(model3_imp)),
  Model4_Coeff = c("-", round(model4_imp_coeffs["gap", "Estimate"], 2), "Yes",
                   round(model4_imp_sigma, 2), nobs(model4_imp))
)

knitr::kable(table4_imputed) %>%
  kableExtra::kable_classic(full_width = FALSE)

#==================================================
# 4.4. Compare the primary estimate against the imputed-sample estimate
#==================================================
primary_path <- here("data", "classic_did_summary.csv")
if (file.exists(primary_path)) {
  primary <- read_csv(primary_path, show_col_types = FALSE)
  cat("\n=== Gap coefficient: primary (complete-case) vs imputed sample ===\n")
  cat(sprintf("Primary (model5, complete-case): %.4f\n", primary$estimate))
  cat(sprintf("Imputed sample (model4):         %.4f\n", model4_imp_coeffs["gap", "Estimate"]))
} else {
  cat("\n(Run 03_classic_did_regression.R first to compare against the primary estimate.)\n")
}

#==================================================
# 4.5. Save outputs
#==================================================
write_csv(est_df_imputed, here("data", "estimation_sample_robustness_check.csv"))
write_csv(table4_imputed, here("output", "tables", "table4_robustness_check.csv"))

cat("\nSaved: data/estimation_sample_imputed.csv, output/tables/table4_sensitivity_imputed.csv\n")