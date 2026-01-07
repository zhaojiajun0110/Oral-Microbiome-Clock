# ==============================================================================
# Script: 03_OMAA_Disease_Incremental_Value.R
# Purpose: Evaluate OMAA incremental value for chronic diseases.
# Output: Delta AUC and Bootstrap P values for the Validation Set.
# ==============================================================================

library(caret)
library(pROC)
library(dplyr)
library(randomForest)
library(openxlsx)
library(tidyverse)

# 1. Variable Definitions ------------------------------------------------------
outcomes <- c(
  "Psoriasis", "Arthritis", "Gout", "Congestive_heart_failure", 
  "Coronary_heart_disease", "Angina", "Heart_disease", "Stroke", 
  "Emphysema", "Liver", "Chronic_bronchitis", "Cancer", "thyroid",
  "diabetes", "hpd"
)

# Baseline covariates
covariates <- c(
  "age", "gender", "race", "educational_level", "marriage_status", 
  "PIR", "drinking_status", "smoking_status", "sleep", "BMI", "met", "tooth count"
)

# Model Predictors
predictors_A <- covariates
predictors_B <- c("aa", covariates) # 'aa' is OMAA

# 2. Training (Discovery Set: NHANES 2009-2010) ---------------------------------
ctrl <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 5, # Optimized repeats for efficiency
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

trained_models <- list()

for (outcome in outcomes) {
  # Prep Discovery Data
  data_train <- merge_09 %>% 
    dplyr::select(aa, all_of(c(outcome, covariates))) %>% 
    drop_na()
  
  data_train[[outcome]] <- factor(ifelse(data_train[[outcome]] == 1, "yes", "no"))
  
  # Train Model A (Baseline)
  set.seed(123)
  m_A <- train(
    as.formula(paste(outcome, "~", paste(predictors_A, collapse = "+"))),
    data = data_train, method = "rf", metric = "ROC", trControl = ctrl,
    ntree = 400, tuneGrid = data.frame(mtry = floor(sqrt(length(predictors_A))))
  )
  
  # Train Model B (Baseline + OMAA)
  set.seed(123)
  m_B <- train(
    as.formula(paste(outcome, "~", paste(predictors_B, collapse = "+"))),
    data = data_train, method = "rf", metric = "ROC", trControl = ctrl,
    ntree = 400, tuneGrid = data.frame(mtry = floor(sqrt(length(predictors_B))))
  )
  
  trained_models[[outcome]] <- list(model_A = m_A, model_B = m_B)
  message("Trained models for: ", outcome)
}

# 3. Validation and P-value Calculation (Validation Set: NHANES 2011-2012) ------
val_summary <- list()

for (outcome in outcomes) {
  # Prep Validation Data
  data_val <- merge_11 %>% 
    dplyr::select(aa, all_of(c(outcome, covariates))) %>% 
    drop_na()
  
  data_val[[outcome]] <- factor(ifelse(data_val[[outcome]] == 1, "yes", "no"))
  
  # Predict Probabilities
  prob_A <- predict(trained_models[[outcome]]$model_A, newdata = data_val, type = "prob")[, "yes"]
  prob_B <- predict(trained_models[[outcome]]$model_B, newdata = data_val, type = "prob")[, "yes"]
  
  # ROC objects
  roc_A <- roc(data_val[[outcome]], prob_A, quiet = TRUE)
  roc_B <- roc(data_val[[outcome]], prob_B, quiet = TRUE)
  
  # Bootstrap Test for P-value (2000 iterations)
  set.seed(123)
  boot_res <- roc.test(roc_A, roc_B, method = "bootstrap", boot.n = 2000, paired = TRUE)
  
  # Store Results: Only Outcome, Delta AUC, and Bootstrap P
  val_summary[[outcome]] <- tibble(
    Outcome     = outcome,
    Delta_AUC   = as.numeric(roc_B$auc - roc_A$auc),
    Bootstrap_P = boot_res$p.value
  )
}

# 4. Export Final Table --------------------------------------------------------
final_table <- bind_rows(val_summary)
write.xlsx(final_table, "results/Validation_AUC_Comparison.xlsx")

message("Validation results exported successfully.")
