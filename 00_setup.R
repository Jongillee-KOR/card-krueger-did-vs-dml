#==================================================
# 0. Setup: Install & Load Packages
#==================================================

required_pkgs <- c("readxl", "dplyr", "tidyr", "lubridate", "readr", "stringr",
                   "here", "ggplot2", "broom", "lmtest", "sandwich",
                   "glmnet", "xgboost", "kableExtra")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) install.packages(new_pkgs, repos = "https://cloud.r-project.org")
invisible(lapply(required_pkgs, library, character.only = TRUE))

dir.create(here("output", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "tables"),  recursive = TRUE, showWarnings = FALSE)