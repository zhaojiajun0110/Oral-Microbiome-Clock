# ==============================================================================
# Script: OMAA_Predictors_Diet_and_Medication.R
# Purpose: Assess the influence of dietary patterns and medication on OMAA.
# Methodology: Random Forest Regression, 10-fold CV, and SHAP value analysis.
# ==============================================================================

library(randomForest)
library(caret)
library(fastshap)
library(openxlsx)
library(tidyverse)

# 1. Data Loading and Preparation ----------------------------------------------

# Datasets should contain SEQN, OMAA (aa), Dietary variables (DR1T...), 
# and binary Medication variables.
df_09 <- read.xlsx("data/contribution_09.xlsx")
df_11 <- read.xlsx("data/contribution_11.xlsx")

# Identify feature sets
diet_vars <- grep("^DR1T", colnames(df_09), value = TRUE)
med_vars  <- setdiff(colnames(df_09), c("SEQN", "aa", diet_vars))

# 2. Dietary Pattern Analysis --------------------------------------------------

set.seed(123)
ctrl <- trainControl(method = "cv", number = 10, savePredictions = "all")

# Train Random Forest for Dietary Patterns
set.seed(123)
rf_diet <- train(
  as.formula(paste("aa ~", paste(diet_vars, collapse = " + "))),
  data = df_09, method = "rf", trControl = ctrl, ntree = 500
)

# Cross-validated Spearman correlation (Discovery)
diet_cv_preds <- rf_diet$pred %>% arrange(rowIndex)
cor_diet_dis <- cor.test(diet_cv_preds$obs, diet_cv_preds$pred, method = "spearman")

# External Validation (Validation set)
diet_val_preds <- predict(rf_diet, newdata = df_11)
cor_diet_val <- cor.test(df_11$aa, diet_val_preds, method = "spearman")

# 3. Medication Use Analysis ---------------------------------------------------

# Train Random Forest for Medication Use
set.seed(123)
rf_med <- train(
  as.formula(paste("aa ~", paste(med_vars, collapse = " + "))),
  data = df_09, method = "rf", trControl = ctrl, ntree = 500
)

# Cross-validated Spearman correlation (Discovery)
med_cv_preds <- rf_med$pred %>% arrange(rowIndex)
cor_med_dis <- cor.test(med_cv_preds$obs, med_cv_preds$pred, method = "spearman")

# External Validation (Validation set)
med_val_preds <- predict(rf_med, newdata = df_11)
cor_med_val <- cor.test(df_11$aa, med_val_preds, method = "spearman")

# 4. Global Feature Importance (SHAP Analysis) ---------------------------------

# Prediction wrapper for fastshap
pred_wrapper <- function(object, newdata) predict(object, newdata = newdata)

# Calculate SHAP values for Medication (Global Importance)
set.seed(123)
shap_med <- explain(
  object = rf_med$finalModel, 
  X = df_09[, med_vars], 
  pred_wrapper = pred_wrapper, 
  nsim = 100
)

shap_summary <- shap_med %>%
  as.data.frame() %>%
  summarise(across(everything(), ~ mean(abs(.)))) %>%
  pivot_longer(cols = everything(), names_to = "Feature", values_to = "MeanAbsSHAP") %>%
  arrange(desc(MeanAbsSHAP))

# 5. Export Results ------------------------------------------------------------

# Summary of model performance
performance_summary <- data.frame(
  Analysis = c("Dietary", "Medication"),
  Discovery_Rho = c(cor_diet_dis$estimate, cor_med_dis$estimate),
  Discovery_P = c(cor_diet_dis$p.value, cor_med_dis$p.value),
  Validation_Rho = c(cor_diet_val$estimate, cor_med_val$estimate),
  Validation_P = c(cor_diet_val$p.value, cor_med_val$p.value)
)

write.xlsx(performance_summary, "results/OMAA_Contribution_Performance.xlsx")
write.xlsx(shap_summary, "results/Medication_SHAP_Importance.xlsx")

message("Contribution analysis complete. Performance metrics and SHAP values exported.")
