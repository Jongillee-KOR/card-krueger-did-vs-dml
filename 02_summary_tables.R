#==================================================
# 2. Summary Tables (Table 1-3)
#==================================================
source(here::here("00_setup.R"))

#==================================================
# 2.1. Load cleaned data
#==================================================
df <- read_csv(here("data", "fastfood_data.csv"), show_col_types = FALSE) %>%
  mutate(state_name = if_else(state == 1, "NJ", "PA"))

#==================================================
# 2.2. Table 3: mean FTE before/after, by state
#==================================================

# Row 1,2: mean and SE by state
row12 <- df %>%
  group_by(state_name) %>%
  summarise(
    mean_before = mean(fte, na.rm = TRUE),
    se_before   = sd(fte, na.rm = TRUE) / sqrt(sum(!is.na(fte))),
    mean_after  = mean(fte2, na.rm = TRUE),
    se_after    = sd(fte2, na.rm = TRUE) / sqrt(sum(!is.na(fte2))),
    .groups = "drop")

# Row 3: change in mean FTE
row3 <- row12 %>%
  mutate(
    mean_diff = mean_after - mean_before,
    se_diff   = sqrt(se_before^2 + se_after^2))

# Row 4: balanced sample
row4 <- df %>%
  filter(!is.na(fte), !is.na(fte2)) %>%
  mutate(diff = fte2 - fte) %>%
  group_by(state_name) %>%
  summarise(
    mean_diff = mean(diff, na.rm = TRUE),
    se_diff   = sd(diff, na.rm = TRUE) / sqrt(sum(!is.na(diff))),
    .groups = "drop")

# Row 5: closed stores set to FTE = 0 (status2 == 3 means store closed)
row5 <- df %>%
  mutate(fte2_adj = if_else(status2 == 3, 0, fte2)) %>%
  filter(!is.na(fte), !is.na(fte2_adj)) %>%
  mutate(diff = fte2_adj - fte) %>%
  group_by(state_name) %>%
  summarise(
    mean_diff = mean(diff, na.rm = TRUE),
    se_diff   = sd(diff, na.rm = TRUE) / sqrt(sum(!is.na(diff))),
    .groups = "drop")

get_row <- function(df, mean_col, se_col) {
  pa <- df %>% filter(state_name == "PA")
  nj <- df %>% filter(state_name == "NJ")
  c(
    PA_Mean = pa[[mean_col]],
    PA_SE = pa[[se_col]],
    NJ_Mean = nj[[mean_col]],
    NJ_SE = nj[[se_col]],
    NJ_PA_Mean = nj[[mean_col]] - pa[[mean_col]],
    NJ_PA_SE = sqrt(nj[[se_col]]^2 + pa[[se_col]]^2))
}

r1 <- get_row(row12, "mean_before", "se_before")
r2 <- get_row(row12, "mean_after", "se_after")
r3 <- get_row(row3, "mean_diff", "se_diff")
r4 <- get_row(row4, "mean_diff", "se_diff")
r5 <- get_row(row5, "mean_diff", "se_diff")

table3 <- data.frame(
  Variable = c(
    "FTE employment before, all observations",
    "FTE employment after, all observations",
    "Change in mean FTE employment",
    "Change in mean FTE employment, balanced sample of stores",
    "Change in mean FTE employment setting FTE at 0 for closed stores"),
  rbind(r1, r2, r3, r4, r5),
  row.names = NULL) %>%
  mutate(across(where(is.numeric), ~round(.x, 2)))

knitr::kable(table3) %>%
  kableExtra::kable_classic(full_width = FALSE)

#==================================================
# 2.3. Table 1: summary statistics
#==================================================
summary_statistics <- df %>%
  group_by(state_name) %>%
  summarise(
    mean_before = mean(fte, na.rm = TRUE),
    mean_after  = mean(fte2, na.rm = TRUE),
    var_before  = var(fte, na.rm = TRUE),
    var_after   = var(fte2, na.rm = TRUE),
    count_before = sum(!is.na(fte)),
    count_after  = sum(!is.na(fte2)),
    .groups = "drop"
  ) %>%
  mutate(
    se_before = sqrt(var_before / count_before),
    se_after  = sqrt(var_after / count_after))

# digits vector length must match number of columns (state_name + 8 numeric cols = 9)
knitr::kable(
  summary_statistics,
  digits = c(NA, 2, 2, 2, 2, 0, 0, 2, 2)
) %>%
  kableExtra::kable_classic(full_width = FALSE)

#==================================================
# 2.4. Table 2: NJ vs PA comparison
#==================================================
PA <- summary_statistics %>% filter(state_name == "PA")
NJ <- summary_statistics %>% filter(state_name == "NJ")

NJ_PA_before    <- NJ$mean_before - PA$mean_before
NJ_PA_before_se <- sqrt(NJ$se_before^2 + PA$se_before^2)

NJ_PA_after    <- NJ$mean_after - PA$mean_after
NJ_PA_after_se <- sqrt(NJ$se_after^2 + PA$se_after^2)

did_mean <- NJ_PA_after - NJ_PA_before
did_se   <- sqrt(NJ_PA_after_se^2 + NJ_PA_before_se^2)

comparison <- data.frame(
  Variable = c(
    "FTE employment before",
    "SE employment before",
    "FTE employment after",
    "SE employment after",
    "Change in mean FTE",
    "SE of change in mean"),
  NJ_PA = c(
    NJ_PA_before,
    NJ_PA_before_se,
    NJ_PA_after,
    NJ_PA_after_se,
    did_mean,
    did_se))

knitr::kable(comparison, digits = 2) %>%
  kableExtra::kable_classic(full_width = FALSE)

#==================================================
# 2.5. Save outputs
#==================================================
write_csv(summary_statistics, here("output", "tables", "table1_summary_statistics.csv"))
write_csv(table3,             here("output", "tables", "table3_replication.csv"))
write_csv(comparison,         here("output", "tables", "table2_nj_pa_comparison.csv"))

cat("Saved: output/tables/table1_summary_statistics.csv, table3_replication.csv, table2_nj_pa_comparison.csv\n")