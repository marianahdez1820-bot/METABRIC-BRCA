#> In this script we obtain and preprocess the metadata.gse96058 and the counts data

library(GEOquery)
library(tidyverse)

# 1.- Download data -----------------------------------------------------------


# 1.1 Download the supplementary file (The actual expression matrix)
# getGEOSuppFiles("GSE96058", baseDir = "D:/GSE96058")

# 1.1.2 Read the specific expression file (SCAN-B typically provides a large .txt or .csv)

raw <- read.csv("D:/GSE96058/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz", 
                row.names = 1, check.names = FALSE)
counts_data.gse96058 <- raw
# 1.2 Download metadata.gse96058

gse <- getGEO("GSE96058", GSEMatrix = TRUE)

# 1.2.2 Asign to object

pheno <- pData(gse[[1]])



# 2.- Preprocess metadata.gse96058 -------------------------------------------------

pheno$characteristics_ch1.3 <- gsub("\\D", "", pheno$characteristics_ch1.3)


metadata.gse96058 <-
  pheno %>%
  mutate(
    tissue = source_name_ch1,
    age = as.numeric(`age at diagnosis:ch1`),
    tumor_size = as.numeric(`tumor size:ch1`),
    lymph_group = `lymph node group:ch1`,
    lymph_status =  `lymph node status:ch1`,
    er_status = as.numeric(`er status:ch1`),
    pgr_status = as.numeric(`pgr status:ch1`),
    her2_status = as.numeric(`her2 status:ch1`),
    ki67_status = as.numeric(`ki67 status:ch1`),
    nhg = as.factor(`nhg:ch1`),
    er_pred_mgc = as.numeric(`er prediction mgc:ch1`),
    # This are predictions made by RNA if mgc its Molecular Gene Classifier which is older than SCN which is Single Sample Classifier (SCAN-B) and if SGC its single gene classifier
    er_pred_sgc = as.numeric(`er prediction sgc:ch1`),
    pgr_pred_mfc = as.numeric(`pgr prediction mgc:ch1`),
    pgr_pred_sgc = as.numeric(`pgr prediction sgc:ch1`),
    her2_pred_mfc = as.numeric(`her2 prediction mgc:ch1`),
    her2_pred_sgc = as.numeric(`her2 prediction sgc:ch1`),
    ki67_pred_mfc = as.numeric(`ki67 prediction mgc:ch1`),
    ki67_pred_sgc = as.numeric(`ki67 prediction sgc:ch1`),
    nhg_pred_mgc = as.numeric(`nhg prediction mgc:ch1`),
    pam50 = as.factor(`pam50 subtype:ch1`),
    os_months = as.numeric(`overall survival days:ch1`) / 30.4166667,
    os_status = as.numeric(`overall survival event:ch1`),
    endocrine_tx = as.numeric(`endocrine treated:ch1`),
    chemo_tx = as.numeric(`chemo treated:ch1`)
    
    
  ) %>%
  dplyr::select(
    -c(
      source_name_ch1,
      characteristics_ch1.2,
      characteristics_ch1.3,
      characteristics_ch1.4,
      characteristics_ch1.5,
      characteristics_ch1.6,
      characteristics_ch1.7,
      characteristics_ch1.8,
      characteristics_ch1.9,
      characteristics_ch1.10,
      characteristics_ch1.11,
      characteristics_ch1.12,
      characteristics_ch1.13,
      characteristics_ch1.14,
      characteristics_ch1.15,
      characteristics_ch1.16,
      characteristics_ch1.17,
      characteristics_ch1.18,
      characteristics_ch1.19,
      characteristics_ch1.20,
      characteristics_ch1.21,
      characteristics_ch1.22,
      characteristics_ch1.23,
      characteristics_ch1.24,
      # Each one of the characteristics_ch1. corresponds to its equivalent in the next lines and both correspond in orther to its characteristic in mutate
      `age at diagnosis:ch1`,
      `tumor size:ch1`,
      `lymph node group:ch1`,
      `lymph node status:ch1`,
      `er status:ch1`,
      `pgr status:ch1`,
      `her2 status:ch1`,
      `ki67 status:ch1`,
      `nhg:ch1`,
      `er prediction mgc:ch1`,
      `er prediction sgc:ch1`,
      `pgr prediction mgc:ch1`,
      `pgr prediction sgc:ch1`,
      `her2 prediction mgc:ch1`,
      `her2 prediction sgc:ch1`,
      `ki67 prediction mgc:ch1`,
      `ki67 prediction sgc:ch1`,
      `nhg prediction mgc:ch1`,
      `pam50 subtype:ch1`,
      `overall survival days:ch1`,
      `overall survival event:ch1`,
      `endocrine treated:ch1`,
      `chemo treated:ch1`
      
    )
  )

counts_data.gse96058[1:5,1:5]


metadata.gse96058_er_pos <-
  metadata.gse96058 %>% 
  filter(er_status == 1) %>% 
  rownames_to_column("id") %>% 
  mutate(EVENT_STAT = os_status,
         EVENT_MON = os_months,
         id = NULL) 

# 3.- Preprocess data -----------------------------------------------------

# 3.1 Match it with metadata.gse96058

# 3.1.1 Identify patients in both sets

common_samples <- intersect(colnames(counts_data.gse96058), metadata.gse96058_er_pos$title)

# 3.1.2 Keep the patients in counts data that also have metadata.gse96058

counts_data.gse96058_erpos <- counts_data.gse96058[, common_samples]

counts_data.gse96058_erpos <- t(counts_data.gse96058_erpos)



# 3.2 Find genes present in both data sets

colnames(counts_data.gse96058_erpos) <- make.names(colnames(counts_data.gse96058_erpos))

common_genes_meta.gse96058 <- intersect(proof_genes, colnames(counts_data.gse96058_erpos))


# 4.2 Object with all the patients and expression of only the genes of interest

counts_data.gse96058_erpos <- counts_data.gse96058_erpos[, colnames(counts_data.gse96058_erpos) %in% common_genes_meta.gse96058]

# 4.3 Stop running if there are less genes in TCGA than on the signature

if(length(common_genes_meta.gse96058) < length(proof_genes)){
  stop(paste("There are missing genes in TCGA relative to the signature, missing ", length(proof_genes) - length(common_genes_meta.gse96058), " gene(s): "),  paste0(proof_genes[!(proof_genes %in% common_genes_meta.gse96058)], sep = ", ")) # Script stops here
}else{
  print("All genes in the signature are on TCGA")
}




counts_data.gse96058_erpos <- counts_data.gse96058_erpos[ , common_genes_meta.gse96058]

proof_genes_pt.gse96058 <- 
  counts_data.gse96058_erpos %>% 
  as.data.frame() %>% 
  rownames_to_column("title") %>% 
  left_join(metadata.gse96058_er_pos, by = "title") %>% 
  mutate(surv_obj = Surv(time = EVENT_MON, event = EVENT_STAT, type = "right")) %>% 
  dplyr::select(all_of(proof_genes),
                EVENT_STAT,
                EVENT_MON,
                title,
                surv_obj) %>% 
  column_to_rownames("title")

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

# 6.- Validation ----------------------------------------------------------


# 6.1 Extract the trained recipe from the workflow

trained_rec <- extract_recipe(final_fit)

# 6.2 "Bake" the RNA-seq data and by that we mean to apply the same steps of the recipe to the new data

proof_genes_pt.gse96058_baked <- bake(trained_rec, new_data = proof_genes_pt.gse96058)

# 6.3 This is the predict phase

gse96058_results <- predict(final_fit, new_data = proof_genes_pt.gse96058_baked, type = "linear_pred") %>%
  bind_cols(proof_genes_pt.gse96058_baked)


# 6.3 To get the p value

validation_test <- coxph(Surv(EVENT_MON, EVENT_STAT) ~ .pred_linear_pred, data = gse96058_results)

summary(validation_test)


library(survminer)

# Create risk groups based on the median of the predictions

gse96058_results <- gse96058_results %>%
  mutate(risk_group = as.factor(ifelse(.pred_linear_pred < median(.pred_linear_pred), "High Risk", "Low Risk")))


# Fit the KM curve

km_fit <- survfit(Surv(EVENT_MON, EVENT_STAT) ~ risk_group, data = gse96058_results)

# Plot

ggsurvplot(km_fit, 
           data = gse96058_results, 
           pval = TRUE, 
           risk.table = TRUE,
           title = "Validation in GSE2034 (Untreated Cohort)",
           palette = c("#E41A1C", "#377EB8"))


gse96058_results <- gse96058_results %>%
  mutate(pred_z = scale(.pred_linear_pred))

# Run Cox again
gse96058_results$risk_group <- relevel(gse96058_results$risk_group, ref = "Low Risk")

summary_gse96058 <- summary(coxph(Surv(EVENT_MON, EVENT_STAT) ~ risk_group, data = gse96058_results))

# Calculate the actual Concordance Index
c_index_results.gse96058 <- concordance(Surv(EVENT_MON, EVENT_STAT) ~ .pred_linear_pred, 
                                        data = gse96058_results)



library(timeROC)

# Area under the curve per time

res_auc <- timeROC(T = gse96058_results$EVENT_MON,
                   delta = gse96058_results$EVENT_STAT,
                   marker = -gse96058_results$pred_z,
                   cause = 1, # The event code
                   times = c(36, 60, 72, 80), # 3, 5, and 6 years
                   iid = TRUE)

# View the AUC values

res_auc_gse96058 <- res_auc$AUC %>% 
  as.data.frame()


# Cox multivariado con clinica

proof_genes_pt_gse96058.cox <- 
  proof_genes_pt.gse96058 %>% 
  as.data.frame() %>% 
  rownames_to_column("title") %>% 
  left_join(metadata.gse96058_er_pos, by = "title", suffix = c("", ".y")) %>%
  dplyr::select(-ends_with(".y")) %>% 
  column_to_rownames("title") %>% 
  mutate(HER2 = her2_pred_sgc, 
         LYMPH = lymph_group,
         PAM50 = pam50,
         AGE = as.numeric(age),
         KI67 = ki67_pred_sgc,
         SCORE = gse96058_results$.pred_linear_pred,
         RISK = gse96058_results$risk_group
  ) %>% 
  dplyr::select(all_of(proof_genes),
                surv_obj,
                LYMPH,
                PAM50,
                AGE,
                HER2,
                KI67,
                SCORE,
                RISK
  ) %>% 
  na.omit()

independent_prog.gse96058 <- coxph(surv_obj ~ PAM50 + KI67 + HER2 + AGE + LYMPH + SCORE, 
                                   data = proof_genes_pt_gse96058.cox) %>% 
  tidy(exponentiate = TRUE, conf.int = TRUE)

