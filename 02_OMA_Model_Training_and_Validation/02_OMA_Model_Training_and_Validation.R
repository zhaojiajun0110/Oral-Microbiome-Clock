# ==============================================================================
# Script: 02_OMA_Clock_Construction.R
# Purpose: Construct the Oral Microbiome Age (OMA) clock using Random Forest.
# Analysis: 10-fold CV in Discovery and validation in two external cohorts.
# ==============================================================================

library(tidyverse)
library(randomForest)
library(caret)
library(openxlsx)

# 1. Data Preparation ----------------------------------------------------------

# Load filtered features and datasets
# Datasets should contain SEQN, age (target), and CLR-transformed genera
discovery_df <- read.xlsx("results/df_09_filtered_features.xlsx")
validation_df_1 <- read.xlsx("results/df_11_filtered_features.xlsx")

# Define predictors (X) and target (y)
X_dis <- discovery_df %>% dplyr::select(-SEQN, -age)
y_dis <- discovery_df$age

# 2. 10-fold Cross-Validation (Discovery Set) ----------------------------------
set.seed(123)
folds <- createFolds(y_dis, k = 10, list = TRUE, returnTrain = FALSE)
cv_preds <- rep(NA, length(y_dis))

for(i in seq_along(folds)){
  test_idx <- folds[[i]]
  train_idx <- setdiff(seq_along(y_dis), test_idx)
  
  rf_tmp <- randomForest(
    x = X_dis[train_idx, ],
    y = y_dis[train_idx],
    ntree = 500,
    mtry = floor(sqrt(ncol(X_dis))),
    importance = TRUE
  )
  cv_preds[test_idx] <- predict(rf_tmp, X_dis[test_idx, ])
}

# Save Discovery CV results
discovery_results <- data.frame(
  SEQN = discovery_df$SEQN,
  Age_Observed = y_dis,
  Age_Predicted = cv_preds
)
write.xlsx(discovery_results, "results/OMA_Discovery_CV_Results.xlsx")

# 3. External Validation 1 (NHANES 2011-2012) ----------------------------------

# Train final model on full Discovery set
set.seed(123)
final_model <- randomForest(
  x = X_dis,
  y = y_dis,
  ntree = 500,
  mtry = floor(sqrt(ncol(X_dis))),
  importance = TRUE
)

# Predict on Validation 1
X_val1 <- validation_df_1 %>% dplyr::select(-SEQN, -age)
val1_preds <- predict(final_model, X_val1)

val1_results <- data.frame(
  SEQN = validation_df_1$SEQN,
  Age_Observed = validation_df_1$age,
  Age_Predicted = val1_preds
)
write.xlsx(val1_results, "results/OMA_Validation_1_Results.xlsx")

# 4. Independent External Validation 2 (Intersection Features) -----------------

ext_val_2_raw <- read.xlsx("data/external_validation_2.xlsx") %>%
  filter(age >= 30)

# Identify common features between Discovery and External Validation 2
dis_features <- colnames(X_dis)
ext2_features <- setdiff(colnames(ext_val_2_raw), c("ID", "age"))
common_features <- intersect(dis_features, ext2_features)

# Retrain model using only common features
X_dis_common <- X_dis %>% dplyr::select(all_of(common_features))
set.seed(123)
model_common <- randomForest(
  x = X_dis_common,
  y = y_dis,
  ntree = 500,
  mtry = floor(sqrt(ncol(X_dis_common)))
)

# Predict on External Validation 2
X_ext2 <- ext_val_2_raw %>% dplyr::select(all_of(common_features))
ext2_preds <- predict(model_common, X_ext2)

ext2_results <- data.frame(
  ID = ext_val_2_raw$ID,
  Age_Observed = ext_val_2_raw$age,
  Age_Predicted = ext2_preds
)
write.xlsx(ext2_results, "results/OMA_External_Validation_2_Results.xlsx")

# 5. Model Summary Output ------------------------------------------------------
cat("OMA Clock Construction Complete.\n")
cat("Discovery MAE:", mean(abs(y_dis - cv_preds)), "\n")
cat("Validation 1 MAE:", mean(abs(val1_results$Age_Observed - val1_results$Age_Predicted)), "\n")
cat("External 2 MAE:", mean(abs(ext2_results$Age_Observed - ext2_results$Age_Predicted)), "\n")
