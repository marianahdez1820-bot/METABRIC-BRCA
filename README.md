# Transcriptomic Signatures Predicting Survival and Recurrence in Estrogen Receptor-Positive Breast Cancer: A Machine Learning Analysis of the METABRIC Cohort

This repository contains the pipeline for preprocessing, feature selection, hyperparameter tuning, training, and external validation, as well as misclassification and enrichment analyses for ER-positive breast cancer data.

## Execution Order & Pipeline Structure

The folders and files are ordered in the exact sequence they should be executed, as objects created in earlier scripts are required for subsequent steps.

---

### 1. METABRIC Preprocessing

#### `1. METABRIC preprocessing/1.1 Data preproccesing.R`
This script handles the initial data pipeline:
* **Steps:** Data loading $\rightarrow$ Duplicate gene management $\rightarrow$ Metadata preprocessing.
* **Key Outputs:** Generates the core objects: `metadata.ER_POS_SURV`, `metadata.ER_POS_REC`, and `counts_data`.

---

### 2. Feature Selection

#### `2. Feature selection/2.1 Boruta.R`
Runs the Boruta feature selection algorithm to identify important transcriptomic predictors.

>  **Important:**
> * **Analysis to be made:** The output depends on the object loaded in **Line 9**. Use `metadata.ER_POS_SURV` for the survival analysis pipeline, or switch to `metadata.ER_POS_REC` for the recurrence pipeline.
> * **Gene Selection:** **Line 43 (Section 5.2)** defines the threshold for how many genes will be analyzed.
> * **Parallelization:** Adjust **Line 67 (Section 6.2)** to allocate the number of CPU threads/cores for processing.
> * **Tentative Genes:** From **Line 97** onward, the script saves the results. Note that **Line 116** forces a final decision on tentative genes.

---

### 3. Signature Preparation

>  **Section summary:**
> 3. Signature Preparation contains the files to output the object used in the regression and the best parameters. The order followed should be either (`3. Prepare_linreg_ER/3.1 Prepare_linreg_ER_rec_global.R`) or (`3. Prepare_linreg_ER/3.2 Prepare_linreg_ER_global_surv.R`) depending on the desired analysis followed by (`3. Prepare_linreg_ER/3.3 Hyperparameter_selection.R`)

#### `3. Prepare_linreg_ER/3.1 Prepare_linreg_ER_rec_global.R`
* **Purpose:** This file prepares the data specifically for the **recurrence analysis**.
* **Function:** It imports and structures the specific gene signature derived from the Boruta feature selection step (`2. Feature selection/2.1 Boruta.R`) and utilizes `metadata.ER_POS_REC` to formulate the object needed for the next steps named `proof_genes_pt`.

#### `3. Prepare_linreg_ER/3.2 Prepare_linreg_ER_global_surv.R`
* **Purpose:** Similar to the previous file, this file prepares the data but contrary to the previous file this one is for the **survival analysis**.
* **Function:** It imports and structures the specific gene signature derived from the Boruta feature selection step (`2. Feature selection/2.1 Boruta.R`) and utilizes `metadata.ER_POS_SURV` to formulate the object needed for the next steps named `proof_genes_pt`.

>  **Important:**
> * **Signature Selection:** **Line 16 (Section 1.2)** For both scripts this line is where the desired signature list object is inserted
> * **`proof_genes_pt`:** This object consists of patients in rows and the genes from the signature in the columns. *EVENT_STAT*, *EVENT_MON* and *surv_obj* also form part of the columns and are the outcome columns. This object is used in most of the posterior files and since it doesnt have an intrinsic identifier that specifies if it comes from survival or recurrence its important to be careful about which signature preparation file is selected. This common object is created with the objective of it being usefull in the downstream analysis without having to do major modifications to the scripts.
> * ""Returning to these scripts:** As stated in the paper not all the genes are found on GSE2034 or GSE96058. As seen on the preparation of these datasets they input an object named `common_genes_meta.gse2034` and `common_genes_meta.gse96058` respectively. The idea is that these objects are passed on to line 16 of either file (`3. Prepare_linreg_ER/3.1 Prepare_linreg_ER_rec_global.R`) or (`3. Prepare_linreg_ER/3.2 Prepare_linreg_ER_global_surv.R`) if genes are missing and to run the pipeleine from here once more now with the signature compatible in all datasets.

#### `3. Prepare_linreg_ER/3.3 Hyperparameter_selection.R`
* **Purpose:** This file consist on the selection of hyperparameters for the Cox regression model.
* **Function:** The input is `proof_genes_pt`. No modifications should be made and the important output is the object `best_params` which is a tibble containing the selected parameters forthe regression.
* **Caveat:** If one wants to observe the different options of best parameters and not just stick with the one chosen because of C-index change the function *select_best()* in line 85 (Section 3.2) to 
*show_best()*


---

## Requirements & Environment
*(Optional: Add a brief list of the core R packages required, like Boruta, tidyverse, etc.)*
