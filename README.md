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

> **Important Pipeline Notes:**
> * **Analysis Mode:** The output depends on the object loaded in **Line 9**. Use `metadata.ER_POS_SURV` for survival analysis, or switch to `metadata.ER_POS_REC` for recurrence.
> * **Gene Selection:** **Line 43 (Section 5.2)** defines the threshold for how many genes will be analyzed.
> * **Parallelization:** Adjust **Line 67 (Section 6.2)** to allocate the number of CPU threads/cores.
> * **Tentative Genes:** **Line 116** forces a final decision on tentative genes before saving results.

---

### 3. Signature Preparation

> **Section Workflow:** This section outputs the final data objects and optimal hyperparameters used for downstream regression. **Run either 3.1 OR 3.2** depending on your desired analysis, followed by **3.3**.

#### `3. Prepare_linreg_ER/3.1 Prepare_linreg_ER_rec_global.R` (Recurrence Path)
* **Purpose:** Prepares data specifically for **recurrence analysis**.
* **Function:** Imports the gene signature from Boruta step 2.1 and utilizes `metadata.ER_POS_REC` to create `proof_genes_pt` and `ml_metadata`.

#### `3. Prepare_linreg_ER/3.2 Prepare_linreg_ER_global_surv.R` (Survival Path)
* **Purpose:** Prepares data specifically for **survival analysis**.
* **Function:** Imports the gene signature from Boruta step 2.1 and utilizes `metadata.ER_POS_SURV` to create `proof_genes_pt` and `ml_metadata`.

> **Critical Object Metadata Caveats:**
> * **Signature Selection:** For both scripts, **Line 16 (Section 1.2)** is where you insert your desired signature list object.
> * **`proof_genes_pt`:** This object contains patients in rows and signature genes in columns, alongside outcome columns (*EVENT_STAT*, *EVENT_MON*, and *surv_obj*). Because it does *not* contain an intrinsic label stating whether it represents survival or recurrence, be highly careful about which preparation file you ran. 
> * **`ml_metadata`:** An outcome-neutral renaming of the filtered patient metadata, allowing downstream scripts to run without manual object name changes.
> * ** External Validation Loop (GSE2034 / GSE96058):** As noted in the paper, some signature genes may be missing in external datasets. During validation preparation, pass the generated `common_genes_meta.gse2034` or `common_genes_meta.gse96058` objects back into **Line 16** of either script 3.1 or 3.2 to re-run the pipeline with a fully compatible dataset-wide signature.

#### `3. Prepare_linreg_ER/3.3 Hyperparameter_selection.R`
* **Purpose:** Selects hyperparameters for the Cox regression model using `proof_genes_pt`.
* **Output:** Generates `best_params` (a tibble containing optimal parameters). 
* **Tip:** To inspect alternative hyperparameter options instead of automatically choosing the top C-index, change `select_best()` to `show_best()` on **Line 85 (Section 3.2)**.

---

### 4. Cox Regression

> **Section Workflow:** Subdivided into the main model (`4.1`) and downstream analyses (`4.2`). This section automatically adapts to whichever file you ran back in Section 3, so ensure your active global objects match your intended analysis (Survival vs. Recurrence).

#### `4. Cox Regression/4.1 Cox_regression_global.R`
* **Purpose:** Handles the METABRIC dataset split, model training, and initial testing.
* **Required Inputs:** `proof_genes_pt`, `ml_metadata`, and `best_params` (optional; if missing, the script recalculates the tuning grid, though this diverges from the published workflow).
* **Key Results & Outputs:** 
  * Main paper findings derive from `summary_cox` and `independent_prog`.
  * Generates downstream analysis objects (`proof_genes_pt.cox`) and image generation inputs (`fit_km`, `plot_roc`, `facet_labels`, `cox_p_metabric`) utilized later in folder `7. Images`.

---

## Requirements & Environment
*(Optional: Add a brief list of the core R packages required, like Boruta, tidyverse, etc.)*
