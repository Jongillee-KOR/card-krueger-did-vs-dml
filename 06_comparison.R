#==================================================
# 6. Comparison: Classic DiD vs DML
#==================================================
source(here::here("00_setup.R"))

#==================================================
# 6.1. Load results
#==================================================
classic_path <- here("data", "classic_did_summary.csv")
dml_path     <- here("data", "dml_summary.csv")

if (!file.exists(classic_path)) stop("Could not find ", classic_path, ". Run 03_classic_did_regression.R first.")
if (!file.exists(dml_path))     stop("Could not find ", dml_path,     ". Run 05_dml.R first.")

classic_did <- read_csv(classic_path, show_col_types = FALSE)
dml         <- read_csv(dml_path,     show_col_types = FALSE)

#==================================================
# 6.2. Build comparison table
#==================================================
comparison_table <- bind_rows(classic_did, dml) %>%
  mutate(
    lower_95 = estimate - 1.96 * std_error,
    upper_95 = estimate + 1.96 * std_error
  )

cat("\n=== Comparison: Classic DiD vs DML ===\n")
knitr::kable(comparison_table, digits = 4) %>%
  kableExtra::kable_classic(full_width = FALSE)
print(comparison_table)

#==================================================
# 6.3. Comparison plot (point estimate + 95% CI)
#==================================================
p_comparison <- ggplot(comparison_table, aes(x = method, y = estimate, color = method)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.2, linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  theme_classic() +
  labs(
    title = "Comparison of Minimum Wage Impact: Classic DiD vs. DML",
    subtitle = "Estimated Effect on Employment Change (Delta FTE)",
    x = "Estimation Method", y = "Coefficient Estimate & 95% CI"
  ) +
  theme(legend.position = "none")

print(p_comparison)
ggsave(here("output", "figures", "classic_did_vs_dml.png"), p_comparison, width = 8, height = 6, dpi = 300)

#==================================================
# 6.4. Heatmap: estimate & SE across methods
#==================================================
heatmap_data <- comparison_table %>%
  select(method, estimate, std_error) %>%
  pivot_longer(cols = c(estimate, std_error), names_to = "metric", values_to = "value")

p_heatmap <- ggplot(heatmap_data, aes(x = method, y = metric, fill = value)) +
  geom_tile(color = "white", linewidth = 1.2) +
  scale_fill_gradient(low = "#EBF3FB", high = "#1D6FA5") +
  geom_text(aes(label = round(value, 4)), color = "black", fontface = "bold", size = 4.5) +
  theme_minimal() +
  labs(title = "Model Comparison Heatmap", x = "Estimation Method", y = "Statistic", fill = "Value") +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 11),
    axis.text  = element_text(size = 10, face = "bold"),
    panel.grid = element_blank()
  )

print(p_heatmap)
ggsave(here("output", "figures", "comparison_heatmap.png"), p_heatmap, width = 8, height = 5, dpi = 300)

#==================================================
# 6.5. Save final outputs
#==================================================
write_csv(comparison_table, here("output", "tables", "final_comparison_table.csv"))
write_csv(heatmap_data,     here("output", "tables", "final_heatmap_data.csv"))

cat("\nSaved: output/tables/final_comparison_table.csv, final_heatmap_data.csv\n")
cat("Saved: output/figures/classic_did_vs_dml.png, comparison_heatmap.png\n")