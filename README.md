# Oral-Microbiome-Clock: Predicting Biological Age and Host Health

This repository contains the R implementation and analytical pipeline for the paper: **"An Oral Microbiome Clock Predicts Biological Age and Host Health"**.

## Overview

We developed and validated a robust **Oral Microbiome Age (OMA) Score**  using machine learning (Random Forest) on 16S rRNA gene profiles from two large NHANES cohorts and one independent dataset consisting of 2,550 samples and . This repository provides the code to replicate the identification of age-associated genera, the construction of the OMA score, and the evaluation of **Oral Microbiome Age Acceleration (OMAA) Score** as a predictor for mortality, frailty, and chronic diseases.

---

## Repository Structure

The code is organized into four main modules following the analytical workflow described in the manuscript:

### 1. `01_GAM_Screening.R`

* **Purpose:** Identification of age-associated oral microbial genera.
* **Method:** Utilizes **Generalized Additive Models (GAMs)** to identify genera with significant linear or non-linear correlations with chronological age, adjusting for gender and race.

### 2. `02_OMA_Model_Training_and_Validation.R`

* **Purpose:** Training the Random Forest age predictor and calculating OMA Score.
* **Content:**
* Model training on the **Discovery set** (NHANES 2009-2010).
* Model validation on the **Validation set** (NHANES 2011-2012) and Independent external validation set (doi: 10.1128/mSystems.00630-19).
* Calculation of **OMAA Score** (the residual of OMA regressed on Chronological Age).
* Feature importance analysis using **SHAP values**.



### 3. `03_OMAA_Disease_Incremental_Value.R`

**Purpose**: Assessing the incremental predictive value of the OMAA score for various chronic diseases.
* **Content**:
* Constructs disease risk prediction models for chronic conditions (e.g., Cardiovascular disease, Cancer) using the **Random Forest** algorithm.
**Model A**: Baseline covariates (including sex, race, BMI, smoking status, etc.) + Chronological Age .
**Model B**: Model A + OMAA Score.

* **Statistical Evaluation**:
* Generates Receiver Operating Characteristic (ROC) curves for both models.
* Calculates the **Area Under the Curve (AUC)** to quantify predictive performance.
* Employs a **Bootstrap test** (2000 iterations) to compare the AUCs of Model A and Model B, providing the  value for the incremental improvement gained by adding OMAA Score.



### 4. `04_OMAA_Predictors_Diet_and_Medication.R`

* **Purpose:** Investigating the influence of external factors (dietary patterns and medication) on the OMAA Score.
* **Content:** Random Forest regression to assess the explanatory power of nutrients, food groups, and medication use on OMAA Score.

---

## Requirements

The analysis was performed in **R**. The following packages are required:

```r
install.packages(c("mgcv", "randomForest", "vegan", "pROC", 
                   "survival", "fastshap", "compositions", "tidyverse"))

```

## Data Availability

Microbiome and clinical data used in this study are publicly available from the **National Health and Nutrition Examination Survey (NHANES)**:

* [NHANES 2009-2010 Cycle](https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?BeginYear=2009)
* [NHANES 2011-2012 Cycle](https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?BeginYear=2011)

