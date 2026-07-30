#==================================================
# 3. Classic DiD Regression (Table 4-5) - PRIMARY analysis
#==================================================
source(here::here("00_setup.R"))

#==================================================
# 3.1. Load estimation sample (complete cases)
#==================================================
data_path <- here("data", "estimation_sample.csv")
if (!file.exists(data_path)) {
  stop("Could not find ", data_path, ". Run 01_data_cleaning.R first.")
}
est_df <- read_csv(data_path, show_col_types = FALSE)
cat(sprintf("Loaded estimation sample. Observations: %d\n", nrow(est_df)))

#==================================================
# 3.2. Estimate regression models
#==================================================
# model1: NJ dummy only
# model2: adds chain + ownership controls
# model3: continuous gap version (no controls)
# model4: gap + chain + ownership controls
# model5: gap + chain + ownership + region controls
model1 <- lm(delta_emp ~ state, data = est_df)
model2 <- lm(delta_emp ~ state + co_owned + chain2 + chain3 + chain4, data = est_df)
model3 <- lm(delta_emp ~ gap, data = est_df)
model4 <- lm(delta_emp ~ gap + co_owned + chain2 + chain3 + chain4, data = est_df)
model5 <- lm(delta_emp ~ gap + co_owned + chain2 + chain3 + chain4 + centralj + northj + pa1, data = est_df)

# Extract coefficients and residual SE
model1_coeffs <- summary(model1)$coefficients
model2_coeffs <- summary(model2)$coefficients
model3_coeffs <- summary(model3)$coefficients
model4_coeffs <- summary(model4)$coefficients

model1_sigma <- summary(model1)$sigma
model2_sigma <- summary(model2)$sigma
model3_sigma <- summary(model3)$sigma
model4_sigma <- summary(model4)$sigma

#==================================================
# 3.3. Nested F-tests: do the added controls matter?
#==================================================
# anova() compares the two nested models directly, testing whether the variables added in model2/model4 are jointly significant
anova_state <- anova(model1, model2)
anova_gap   <- anova(model3, model4)

p_controls_state <- anova_state$`Pr(>F)`[2]
p_controls_gap   <- anova_gap$`Pr(>F)`[2]

cat("\n=== Nested F-test: does adding chain+ownership controls matter? ===\n")
cat(sprintf("State model (model1 vs model2): p = %.4f\n", p_controls_state))
cat(sprintf("Gap model   (model3 vs model4): p = %.4f\n", p_controls_gap))

#==================================================
# 3.4. Table 4: New Jersey Dummy specification
#==================================================
table4_calculated <- data.frame(
  Independent_Var = c(
    "New Jersey Dummy",
    "Initial Wage Gap",
    "Controls for Chain and Ownership",
    "Controls for Region",
    "Standard Error of Regression",
    "P-value for Controls (nested F-test)",
    "Number of Stores in Sample"
  ),
  Model1_Coeff = c(round(model1_coeffs["state", "Estimate"], 2), "-", "No", "No",
                   round(model1_sigma, 2), "-", nobs(model1)),
  Model1_SE = c(round(model1_coeffs["state", "Std. Error"], 2), "-", "No", "No",
                round(model1_sigma, 2), "-", nobs(model1)),
  Model2_Coeff = c(round(model2_coeffs["state", "Estimate"], 2), "-", "Yes", "No",
                   round(model2_sigma, 2), round(p_controls_state, 4), nobs(model2)),
  Model2_SE = c(round(model2_coeffs["state", "Std. Error"], 2), "-", "Yes", "No",
                round(model2_sigma, 2), round(p_controls_state, 4), nobs(model2)),
  Model3_Coeff = c("-", round(model3_coeffs["gap", "Estimate"], 2), "No", "No",
                   round(model3_sigma, 2), "-", nobs(model3)),
  Model3_SE = c("-", round(model3_coeffs["gap", "Std. Error"], 2), "No", "No",
                round(model3_sigma, 2), "-", nobs(model3)),
  Model4_Coeff = c("-", round(model4_coeffs["gap", "Estimate"], 2), "Yes", "No",
                   round(model4_sigma, 2), round(p_controls_gap, 4), nobs(model4)),
  Model4_SE = c("-", round(model4_coeffs["gap", "Std. Error"], 2), "Yes", "No",
                round(model4_sigma, 2), round(p_controls_gap, 4), nobs(model4))
)

knitr::kable(table4_calculated) %>%
  kableExtra::kable_classic(full_width = FALSE)

#==================================================
# 3.5. Table 5: Region-controlled specification (model5)
#==================================================
model5_coeffs <- summary(model5)$coefficients
model5_sigma  <- summary(model5)$sigma

table5_regression <- data.frame(
  Independent_Var = c(
    "New Jersey Dummy",
    "Initial Wage Gap",
    "Controls for Chain and Ownership",
    "Controls for Region",
    "Standard Error of Regression",
    "Number of Stores in Sample"
  ),
  Model1_Coeff = c(round(model1_coeffs["state", "Estimate"], 2), "-", "No", "No", round(model1_sigma, 2), nobs(model1)),
  Model4_Coeff = c("-", round(model4_coeffs["gap", "Estimate"], 2), "Yes", "No", round(model4_sigma, 2), nobs(model4)),
  Model5_Coeff = c("-", round(model5_coeffs["gap", "Estimate"], 2), "Yes", "Yes", round(model5_sigma, 2), nobs(model5)),
  Model5_SE    = c("-", round(model5_coeffs["gap", "Std. Error"], 2), "Yes", "Yes", round(model5_sigma, 2), nobs(model5))
)

knitr::kable(table5_regression) %>%
  kableExtra::kable_classic(full_width = FALSE)

#==================================================
# 3.6. Save outputs
#==================================================
write_csv(table4_calculated, here("output", "tables", "table4_classic_did.csv"))
write_csv(table5_regression, here("output", "tables", "table5_region_controls.csv"))

classic_did_summary <- tibble(
  method    = "Classic DiD (OLS, model5)",
  estimate  = model5_coeffs["gap", "Estimate"],
  std_error = model5_coeffs["gap", "Std. Error"]
)
write_csv(classic_did_summary, here("data", "classic_did_summary.csv"))

cat("\nSaved: output/tables/table4_classic_did.csv, table5_region_controls.csv, data/classic_did_summary.csv\n")