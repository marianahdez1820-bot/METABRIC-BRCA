library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)



# 1.- Loading data --------------------------------------------------------



# 1.1 Query con base a RNA-Seq y a 3 casos y 3 controles

# tcga_rna <- GDCquery("TCGA-BRCA",
#                      data.category = "Transcriptome Profiling",
#                      access = "open",
#                      experimental.strategy = "RNA-Seq",
#                      workflow.type = "STAR - Counts"
# )

# GDCdownload(tcga_rna, method = "api", files.per.chunk = 5)

# 1.2 Prepare data for usage

#tcga_brca_data <- GDCprepare(tcga_rna, directory = "D:/tcga/GDCdata")

# 1.3 Count matrix

# brca_matrix <- assay(tcga_brca_data, "fpkm_unstrand")

# brca_matrix %>%
#   write.csv("fpkm_unstrand.csv")

brca_matrix <- read.csv("fpkm_unstrand.csv") %>% 
  column_to_rownames("X")


# 1.4 Convert to data frmae

brca_data <- brca_matrix %>% 
  as.data.frame()


colnames(brca_data) <- gsub("\\.", "-", colnames(brca_data))
# 1.5 Extract sample type from TCGA barcode

# 1.5.2 Select the 14th - 16th value which correspond to sample type codes https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes

sample_type_full <- substr(colnames(brca_data), 14, 16)

# 1.5.3 Maintain only primary tumor samples so as to avoid duplicates

brca_data <- brca_data[, sample_type_full == "01A"]
brca_data_find <- t(brca_data) %>% 
  as.data.frame() %>% 
  rownames_to_column("PATIENT_ID")

brca_data2 <- brca_data

colnames(brca_data2) <- substr(colnames(brca_data2), 1, 15)

# We can see that there are 5 patients with 2 samples of the same tumor

names(brca_data)[substr(colnames(brca_data), 1, 15) %in% names(brca_data2)[duplicated(names(brca_data2))]]

# So we keep only the patient sample with highest variance

brca_data2 <- brca_data

colnames(brca_data2) <- substr(colnames(brca_data2), 1, 15)


sample_variance <- apply(brca_data2, 2, var, na.rm = TRUE)

df <- data.frame(
  sample = colnames(brca_data2),
  variance = sample_variance
)

selected_samples <- df %>%
  group_by(sample) %>%
  slice_max(order_by = variance, n = 1, with_ties = FALSE) %>% 
  pull(sample)

brca_data2 <- brca_data2[, selected_samples]

# 2.- Metadata ------------------------------------------------------------
library(UCSCXenaTools)


# 1. Define the destination directory on your D drive
my_dir <- "D:/tcga/GDCdata/Metadata"
if (!dir.exists(my_dir)) dir.create(my_dir, recursive = TRUE)

# 2. Generate and Query
# We use the specific TCGA BRCA clinical matrix from the public hub
data_query <- XenaGenerate(subset = XenaDatasets == "TCGA.BRCA.sampleMap/BRCA_clinicalMatrix") %>% 
  XenaQuery()

# 3. Download with a specific destination
# We add 'destdir' to ensure it goes to your D drive instead of a temp folder
xe_download <- XenaDownload(data_query, destdir = my_dir)

# 4. Prepare (Load) the data
# Instead of passing a string, we pass the actual download object 
# This prevents the "res not found" error
brca_clinical <- XenaPrepare(xe_download)



brca_clinical$days_to_last_followup



# 1. Create the Recurrence variables
refined_data <- brca_clinical %>%
  mutate(
    # Create the binary EVENT (1 = Yes, 0 = No/Censored)
    event = ifelse(OS_event_nature2012 == "1", 1, 0),
    
    # Calculate the TIME in months
    # Logic: If they recurred, use the 'days_to_new_tumor' column.
    # If they didn't recur, use 'days_to_last_followup' (this is called Censoring).
    event_mon = OS_Time_nature2012,
    # Convert to Months for the - gene model
    event_mon = event_mon / 30.44
  ) %>%
  # 2. Filter for ER status
  filter(
    ER_Status_nature2012 == "Positive", # Only ER+
    event_mon > 0               # Remove invalid/negative entries
  ) 

# Preview the new columns
head(refined_data)


# Count object that corresponds to the metadata patients
brca_data2 <- brca_data2[, colnames(brca_data2) %in% refined_data$sampleID]

# 3.- Deleting duplicates and asigning ensembl as rownames ----------------


# 3.1 Deleting the version of ensembl and keeping only the full name

brca_data2$ensembl <- gsub("\\..*", "", rownames(brca_data2))

 
# 3.2 Variance

numeric_data <- brca_data2 %>%
  dplyr::select(where(is.numeric))

brca_data2$variance <- apply(numeric_data, 1, var)

# 3.2.2 Only mantain the version of the gene duplicate with higher variance

brca_data.unique <- brca_data2 %>% # Initial data
  group_by(ensembl) %>% # Group by ensembl
  slice_max(order_by = variance, n = 1, with_ties = FALSE) %>% # Order by variance and keep the highest
  ungroup() %>% 
  column_to_rownames("ensembl") %>% # Asign ensembl as rownames
  dplyr::select( - variance) # Delete variance column
  
# 3.3 Create object with smbol an genetype

gene_info <- AnnotationDbi::select(org.Hs.eg.db,
                    keys = rownames(brca_data.unique),
                    columns = c("SYMBOL", "GENETYPE"),
                    keytype = "ENSEMBL")

# 3.3.2 Create object with only protein coding genes

protein_coding <- gene_info %>%
  filter(GENETYPE == "protein-coding")


# 3.4 Keep counts of only protein coding genes

brca_data.unique <-  brca_data.unique[
  rownames(brca_data.unique) %in% protein_coding$ENSEMBL, 
]

# 3.5 Asign symbol as a column

brca_data.unique$symbol <- mapIds(
  org.Hs.eg.db,
  keys = rownames(brca_data.unique), # A donde va a buscar
  column = "SYMBOL",  # Nueva columna con ese formato
  keytype = "ENSEMBL",  # Que formato va a buscar en keys
  multiVals = "first") # Que hacer si hay varios del mismo Key




# 4.- Asigning signature genes --------------------------------------------




# 4.2 Object with all the patients and expression of only the genes of interest

late_genes.patients_tcga <- brca_data.unique[brca_data.unique$symbol %in% late_death.genes, ]


# 4.1 List of genes described in the differential expression as being of prognosis for late death

common_genes_meta.tcga <- intersect(boruta_signature, late_genes.patients_tcga$symbol)



late_genes.patients_tcga <- late_genes.patients_tcga %>% 
  mutate(var = matrixStats::rowVars(as.matrix(dplyr::select(., -symbol)))) %>% 
  dplyr::group_by(symbol) %>% 
  slice_max(order_by = var, n = 1, with_ties = FALSE) %>%
  mutate(var = NULL) %>% 
  ungroup() %>% 
  column_to_rownames("symbol")





microarray_data$variance <- microarray_data.var

# 1.2.4 Only maintain the version of the gene duplicate with higher variance

microarray_data.unique <- microarray_data %>% # Initial data
  group_by(Hugo_Symbol) %>% # Group by HUGO symbol
  slice_max(order_by = variance, n = 1) %>% # Order by variance and keep the highest
  ungroup()


# 1. Transpose first so samples are rows
tcga_transposed <- t(late_genes.patients_tcga)

# 2. Log2 Transform (This brings 0-100,000 down to a 0-16 range, like Microarray)
# We add 1 to avoid log(0)
tcga_log <- log2(tcga_transposed + 1)

# 3. Scale (Only if your original 'final_fit' recipe used scaling/normalization)
late_genes.patients_tcga <- tcga_log %>% as.data.frame()


late_genes.patients_tcga <- late_genes.patients_tcga %>%
  as.data.frame() %>% 
  dplyr::select(where(~ !all(is.na(.))))





# 4.4 Check that the patients are in the same order

refined_data <- refined_data[refined_data$sampleID %in% rownames(late_genes.patients_tcga),]


all(rownames(late_genes.patients_tcga) == refined_data$sampleID)

# 4.5 Add a column of EVENT as a binary term for it to be the outcome

late_genes.patients_tcga <- 
  late_genes.patients_tcga %>% 
  rownames_to_column("sampleID") %>% 
  left_join(refined_data, by = "sampleID") %>% 
  column_to_rownames("sampleID") %>% 
  dplyr::select(all_of(late_death.genes), 
         event,
         event_mon) %>% 
  as.data.frame() %>% 
  mutate(EVENT = event, 
         EVENT = as.numeric(event ),
         EVENT_MON = event_mon,
         EVENT_MON = as.numeric(event_mon)
         ) %>%  
  mutate(surv_obj =  Surv(
    time  = EVENT_MON,
    event = EVENT,
    type  = "right"
  )) %>% 
  dplyr::select(- event,
                - event_mon)
late_genes.patients_tcga <- late_genes.patients_tcga %>% 
  filter(EVENT_MON > 0)



###############################################################################
###############################################################################
###############################################################################
###############################################################################
# ==========================================================================
# 8. EXTERNAL VALIDATION ON TCGA
# ==========================================================================


# This is the "testing" phase
tcga_results <- predict(final_fit, new_data = late_genes.patients_tcga, type = "linear_pred") %>%
  bind_cols(late_genes.patients_tcga)


# This is where you get your "P-value" for the test
validation_test <- coxph(surv_obj ~ .pred_linear_pred, data = tcga_results)

summary(validation_test)


library(survminer)

# Create risk groups based on the median of your predictions



tcga_results <- tcga_results %>%
  mutate(risk_group = as.factor(ifelse(.pred_linear_pred < median(.pred_linear_pred), "High Risk", "Low Risk")))

tcga_results$risk_group <- relevel(tcga_results$risk_group, ref = "Low Risk")

# Fit the KM curve
km_fit <- survfit(Surv(EVENT_MON, EVENT) ~ risk_group, data = tcga_results)

# Plot
ggsurvplot(km_fit, 
           data = tcga_results, 
           pval = TRUE, 
           risk.table = TRUE,
           title = "Validation in tcga (Untreated Cohort)",
           palette = c("#E41A1C", "#377EB8"))


tcga_results <- tcga_results %>%
  mutate(pred_z = scale(.pred_linear_pred))

# Run Cox again
summary_cox_tcga <- summary(coxph(Surv(EVENT_MON, EVENT) ~ risk_group, data = tcga_results))

# Calculate the actual Concordance Index
c_index_results.tcga <- concordance(Surv(EVENT_MON, EVENT) ~ .pred_linear_pred, 
                               data = gse96058_results)



library(timeROC)

# Assuming tcga_results has: 
# EVENT_MON (time), EVENT (event), and pred_z (your score)

res_auc <- timeROC(T = tcga_results$EVENT_MON,
                   delta = tcga_results$EVENT,
                   marker = -tcga_results$pred_z,
                   cause = 1, # The event code
                   times = c(36, 60, 120), # 3, 5, and 10 years
                   iid = TRUE)

# View the AUC values
res_auc_tcga <- res_auc$AUC %>% 
  as.data.frame()

# Cox multivariado con clinica
late_genes.patients_tcga.cox <-
  late_genes.patients_tcga %>%
  rownames_to_column("sampleID") %>%
  left_join(refined_data, by = "sampleID") %>%
  as.data.frame() %>%
  column_to_rownames("sampleID") %>%
  mutate(
    HER2 = HER2_Final_Status_nature2012,
    LYMPH = as.numeric(lymph_node_examined_count),
    PAM50 = PAM50Call_RNAseq,
    AGE = as.numeric(Age_at_Initial_Pathologic_Diagnosis_nature2012),
    SCORE = tcga_results$.pred_linear_pred
  ) %>%
  dplyr::select(all_of(late_death.genes),
                surv_obj,
                LYMPH,
                -PAM50,
                AGE,
                -HER2,
                SCORE) %>%
  na.omit()

independent_prog.tcga <- coxph(surv_obj ~ SCORE + AGE + LYMPH, data = late_genes.patients_tcga.cox) %>% 
  tidy(exponentiate = TRUE, conf.int = TRUE)


#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################

text_validation.tcga <- paste("Esta firma en TCGA consiguio un HR de ",
      round(summary_cox_tcga$coefficients[2], 3),
      " (IC 95% de ",
      round(summary_cox_tcga$conf.int[3], 3),
      " - ",
      round(summary_cox_tcga$conf.int[4], 3),
      ", pval de ",
      summary_cox_tcga$coefficients[5],
      ", C score de ",
      round(c_index_results.tcga$concordance, 2),
      ", área bajo la curva a los 3 años de ",
      round(res_auc_tcga[1,1], 3),
      ", a los 5 años de ",
      round(res_auc_tcga[2,1], 3),
      ", y a los 10 años de",
      round(res_auc_tcga[3,1], 3)
)

cat(text_signature,
    text_validation.gse96058,
    text_validation.tcga,
    sep = ". ")





lymph_rows.tcga <- independent_prog.tcga[grepl("LYMPH", independent_prog.tcga$term),]

best_lymph.tcga <- lymph_rows.tcga[which.min(lymph_rows.tcga$p.value),]

significance.tcga <- if((independent_prog.tcga[independent_prog.tcga$term == "SCORE",]$p.value < 0.05) == TRUE){
  "se mantuvo como un predictor de la supervivencia independiente significativo "
}else{
  "no se mantuvo como un predictor de la supervivencia independiente significativo " 
}

text_independent_prog.tcga <- paste0(
  "Con esta base de datos, al realizar un cox multivariado junto a edad y ganglios linfaticos, la firma ",
  significance.tcga,
  "obteniendo un HR de ",
  round(independent_prog.tcga$estimate[independent_prog.tcga$term == "SCORE"], 3),
  " (IC 95% de ",
  round(independent_prog.tcga$conf.low[independent_prog.tcga$term == "SCORE"], 2),
  " - ",
  round(independent_prog.tcga$conf.high[independent_prog.tcga$term == "SCORE"], 2),
  " pvalue de ",
  independent_prog.tcga$p.value[independent_prog.tcga$term == "SCORE"],
  ")",
  if((independent_prog.tcga[independent_prog.tcga$term == "SCORE",]$p.value < best_lymph.tcga$p.value) == TRUE){
    paste0(" superando a los ganglios linfaticos como predictor (HR de ",
           round(best_lymph.tcga$estimate, 3),
           " pval de ",
           best_lymph.tcga$p.value,
           ")")
  }else{
    paste0(" sin lograr superar a los ganglios linfaticos como predictor (HR de ",
           round(best_lymph.tcga$estimate, 3),
           " pval de ",
           best_lymph.tcga$p.value,
           ")")
  }
)

text_tcga <- paste(text_validation.tcga, text_independent_prog.tcga, sep = ". ")

cat(text_metabric, text_gse96058, text_tcga, sep = "\n")
