#==================================================
# 1. Data Cleaning
#==================================================
source(here::here("00_setup.R"))

#==================================================
# 1.1. Load codebook and raw data
#==================================================
codebook <- read_lines(here("data", "codebook"))
dataset  <- read_table(here("data", "public.dat"), col_names = FALSE)

# The codebook rows 8-59 are the hold variable names so I dropped the rows 5,6,13,14,32,33 which are separators
variable_names <- codebook %>%
  .[8:59] %>%
  .[-c(5, 6, 13, 14, 32, 33)] %>%
  str_sub(1, 13) %>%
  str_squish() %>%
  str_to_lower()

dataset <- dataset %>%
  select(-X47) %>%                              # X47 is N/A
  `colnames<-`(., variable_names) %>%
  mutate(across(everything(), as.numeric)) %>%
  mutate(sheet = as.character(sheet))

#==================================================
# 1.2. Create FTE employment and state variables
#==================================================
df <- dataset %>%
  mutate(
    fte  = empft  + nmgrs  + 0.5 * emppt,         # FTE before minimum wage change
    fte2 = empft2 + nmgrs2 + 0.5 * emppt2,        # FTE after minimum wage change
    state_name = if_else(state == 1, "NJ", "PA")
  )

#==================================================
# 1.3. Build the estimation sample 
#==================================================
# keep only stores with all 4 required variables (fte, fte2, wage_st, wage_st2)
est_df <- df %>%
  ungroup() %>%
  filter(!is.na(fte), !is.na(fte2), !is.na(wage_st), !is.na(wage_st2)) %>%
  mutate(delta_emp = fte2 - fte) %>%

    # gap: For New Jersey stores falling below the new $5.05 minimum
    # the wage shortfall ratio is set to 0 for Pennsylvania and high-wage NJ stores.
  
  mutate(gap = if_else(state == 1 & wage_st <= 5.05, (5.05 - wage_st) / wage_st, 0)) %>%
  mutate(
    chain1 = if_else(chain == 1, 1, 0),  # Dummy for chain 1 (Burger King)
    chain2 = if_else(chain == 2, 1, 0),  # Dummy for chain 2 (KFC)
    chain3 = if_else(chain == 3, 1, 0),  # Dummy for chain 3 (Roy Rogers)
    chain4 = if_else(chain == 4, 1, 0)   # Dummy for chain 4 (Wendy's)
  )

cat(sprintf("Full dataset: %d rows. Estimation sample (complete cases): %d rows.\n",
            nrow(df), nrow(est_df)))

#==================================================
# 1.4. Save outputs
#==================================================
write_csv(df,     here("data", "fastfood_data.csv"))
write_csv(est_df, here("data", "estimation_sample.csv"))

cat("Saved: data/fastfood_data.csv, data/estimation_sample.csv\n")