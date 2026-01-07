# ==============================================================================
# Script: GAM_Screening.R
# Purpose: Screen age-associated oral microbial genera using GAM.
# Strategy: Discovery (FDR < 0.05) and Validation (P < 0.05)
# ==============================================================================

library(mgcv)
library(dplyr)
library(openxlsx)
library(tidyverse)

# 1. Load Data -----------------------------------------------------------------
# Ensure 'data/' directory contains dat_09.xlsx and dat_11.xlsx.
df_09 <- read.xlsx("data/dat_09.xlsx") # Discovery: NHANES 2009-2010
df_11 <- read.xlsx("data/dat_11.xlsx") # Validation: NHANES 2011-2012

# Define genera and covariates for the model
genus_cols <- setdiff(names(df_09), c("SEQN", "RIDAGEYR", "gender", "race"))
covar_cols <- c("gender", "race")

# 2. GAM Screening Function ----------------------------------------------------
run_gam_screening <- function(data, genera, covars) {
  results <- list()
  for(g in genera) {
    tryCatch({
      # Escape special characters in genus names
      g_esc <- if(grepl("[-/ ]", g)) sprintf("`%s`", g) else g
      
      # Build formula: Genus ~ s(Age) + Covariates
      form  <- as.formula(sprintf("%s ~ s(RIDAGEYR) + %s", 
                                  g_esc, paste(covars, collapse = " + ")))
      
      # Fit GAM with Restricted Maximum Likelihood (REML)
      fit   <- gam(form, data = data, method = "REML")
      summ  <- summary(fit)
      
      results[[g]] <- data.frame(
        Genus = g, 
        P_val = summ$s.table["s(RIDAGEYR)", "p-value"], 
        EDF   = summ$s.table["s(RIDAGEYR)", "edf"],
        stringsAsFactors = FALSE
      )
    }, error = function(e) message("Error in ", g, ": ", e$message))
  }
  return(bind_rows(results))
}

# 3. Execution & Dual-Cohort Filtering -----------------------------------------

# Discovery Phase: Identify age-associated genera
res_09 <- run_gam_screening(df_09, genus_cols, covar_cols) %>%
  mutate(FDR_09 = p.adjust(P_val, method = "BH")) %>%
  rename(p_09 = P_val, EDF_09 = EDF)

# Validation Phase: Independent verification
res_11 <- run_gam_screening(df_11, genus_cols, covar_cols) %>%
  rename(p_11 = P_val, EDF_11 = EDF)

# Merge and identify stable age-associated genera
sig_genera_list <- res_09 %>%
  inner_join(res_11, by = "Genus") %>%
  filter(FDR_09 < 0.05 & p_11 < 0.05) %>%
  arrange(p_09)

# 4. Export Filtered Feature Matrix --------------------------------------------
# Extract SEQN and selected genera for downstream OMA clock training
final_features <- sig_genera_list$Genus

df_09_sig <- df_09 %>% select(SEQN, RIDAGEYR, all_of(final_features))
df_11_sig <- df_11 %>% select(SEQN, RIDAGEYR, all_of(final_features))

write.xlsx(sig_genera_list, "results/Significant_Genera_Stats.xlsx")
write.xlsx(df_09_sig, "results/df_09_filtered_features.xlsx")
write.xlsx(df_11_sig, "results/df_11_filtered_features.xlsx")

message("GAM screening complete. ", length(final_features), " genera identified.")
