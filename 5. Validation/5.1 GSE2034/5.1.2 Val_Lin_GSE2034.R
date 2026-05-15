
# 6.- Validation ----------------------------------------------------------

# 6.1 Extract the trained recipe from the workflow


# This is the predict phase

results.gse2034 <- predict(final_fit, new_data = proof_genes_pt_gse2034, type = "linear_pred") 

proof_genes_pt_gse2034$.pred_linear_pred <- results.gse2034$.pred_linear_pred

proof_genes_pt_gse2034$id <- metadata.gse.2034_er_pos$file_name

# 6.3 To get the p value

validation_test <- coxph(Surv(EVENT_MON,  EVENT_STAT) ~ .pred_linear_pred, data = proof_genes_pt_gse2034)

summary(validation_test)

# 6.8 Calculate the actual Concordance Index

c_index_results.2034 <- concordance(Surv(EVENT_MON, EVENT_STAT) ~ .pred_linear_pred, 
                                    data = proof_genes_pt_gse2034)

library(survminer)

# Create risk groups based on the median of the predictions

proof_genes_pt_gse2034 <- proof_genes_pt_gse2034 %>%
  mutate(risk_group = as.factor(ifelse(.pred_linear_pred < median(.pred_linear_pred), "High Risk", "Low Risk")))

# Fit the KM curve

km_fit <- survfit(Surv(EVENT_MON,  EVENT_STAT) ~ risk_group, data = proof_genes_pt_gse2034)

# Plot

ggsurvplot(km_fit, 
           data = proof_genes_pt_gse2034, 
           pval = TRUE, 
           risk.table = TRUE,
           title = "Validation in GSE2034 (Untreated Cohort)",
           font.title = 30,
           legend = "bottom",
           font.legend = 22,
           legend.title = "Risk group",
           font.legend.title = 20,
           legend.labs = c("High risk", "Low risk"),
           font.legend.labs = 18,
           xlab = "Time (months)",
           
           xlim = c(0, 180),         # Zoom in
           break.time.by = 50,      # X axis breaks
           ggtheme = theme_minimal(), # ggplot2 theme
           
           linewidth = 3, 
           palette = c("#E41A1C", "#377EB8"),
           )


proof_genes_pt_gse2034 <- proof_genes_pt_gse2034 %>%
  mutate(pred_z = scale(.pred_linear_pred))

# Run Cox again
proof_genes_pt_gse2034$risk_group <- relevel(proof_genes_pt_gse2034$risk_group, ref = "Low Risk")


summary_gse2034 <- summary(coxph(Surv(EVENT_MON,  EVENT_STAT) ~ risk_group, data = proof_genes_pt_gse2034))


library(timeROC)

# Area under the curve per time

res_auc <- timeROC(T = proof_genes_pt_gse2034$EVENT_MON,
                   delta = proof_genes_pt_gse2034$ EVENT_STAT,
                   marker = -proof_genes_pt_gse2034$pred_z,
                   cause = 1, # The event code
                   times = c(36, 60, 120), # 3, 5, and 10 years
                   iid = TRUE)

# 6.9.2 View the AUC values

res_auc_gse2034 <- res_auc$AUC %>% 
  as.data.frame()

# View the AUC values

print(res_auc$AUC)


# Multivariate regression cox with clinical data

proof_genes_pt.gse2034.cox <- 
  proof_genes_pt_gse2034 %>% 
  as.data.frame() %>% 
  rownames_to_column("file_name") %>% 
  left_join(metadata.gse.2034_er_pos, by = "file_name", suffix = c("", ".y")) %>%
  dplyr::select(-ends_with(".y")) %>% 
  column_to_rownames("file_name") %>% 
  mutate(SCORE = proof_genes_pt_gse2034$.pred_linear_pred
  ) %>% 
  dplyr::select(all_of(proof_genes),
                surv_obj,
                SCORE
  ) %>% 
  na.omit()

cox_model.gse2034 <- coxph(surv_obj ~ SCORE, 
                           data = proof_genes_pt.gse2034.cox) 

summary(cox_model.gse2034)

independent_prog.gse2034 <- cox_model.gse2034 %>% 
  tidy(exponentiate = TRUE, conf.int = TRUE)



num_param_compare <- c(9:21)

cat(paste0("The signature got a C-score of ", round(c_index_results.2034$concordance, 2)),
    paste0("an HR of ", round(summary_gse2034$coefficients[2], 2), " (CI 95% of ", round(summary_gse2034$conf.int[3], 2), " - ", round(summary_gse2034$conf.int[4], 2), " pval ", summary_gse2034$coefficients[5], ")"),
    paste0("AUC at 3 years of ", round(res_auc_gse2034[1,], 2), " at 5 years of ", round(res_auc_gse2034[2,], 2), " and at 6 years of "),
    paste0("As an independence factor it has an HR of ", round(independent_prog.gse2034$estimate[independent_prog.gse2034$term == "SCORE"], 2), " (CI 95% of ", round(independent_prog.gse2034$conf.low[independent_prog.gse2034$term == "SCORE"], 2), " - ", round(independent_prog.gse2034$conf.high[independent_prog.gse2034$term == "SCORE"], 2), " pval of ", independent_prog.gse2034$p.value[independent_prog.gse2034$term == "SCORE"], ")"),
    sep = ". "
)





# 4.5 Forest plot ignoring values that tend to infinite

independent_prog.gse2034 %>%
  filter(estimate > 0.0001,
         conf.high < 100) %>%
  mutate(
    term = reorder(term, p.value),
    significant = p.value < 0.05
  ) %>%
  ggplot(aes(x = estimate, y = term, color = significant)) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.5, linewidth = 1.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  theme_minimal()






# 6.- Outlier analysis ----------------------------------------------------


# 6.1 Utilize the fitted object to make predictions based on time

# 6.1.1 Change to numeric 

proof_genes_pt_gse2034 <- 
  proof_genes_pt_gse2034 %>% 
  mutate(EVENT_STAT = as.numeric(EVENT_STAT)) 

# 6.1.2 Augment on different time points

model_diagnostics <- augment(
  final_fit, 
  new_data = proof_genes_pt_gse2034, 
  eval_time = c(36, 60, 120) # Times in months (3, 5, 10 years)
)

# 6.2 Create list to then add the results

list <- list()


#> 6.3 For loop that at each desired time point calculates the bias scores, extracts patients with high bias scores
#> observe metadata of outlier patients and plot the distribution of the prediction with the actual event time

for (i in c(36, 60, 120)) {
  
  eval_time <- i
  
  # 6.3.1 Unnest the predictions and find the biggest outliers on a set point and event
  
  outliers <- model_diagnostics %>%
    dplyr::select(id, EVENT_MON, EVENT_STAT, .pred) %>%
    unnest(.pred) %>%
    filter(.eval_time == eval_time) %>% 
    arrange(desc(.pred_survival)) # Siunce our signature as it goes up, the predicted mortality goes down we see which patients died early who were predicted to die late or survive
  
  
  # 6.3.2 Object with metadata and score characteristics
  
  outlier_summary <- outliers %>%
    inner_join((metadata.gse.2034_er_pos %>% rename(id = "file_name_2", file_name = "id")), by = "id", suffix = c("", ".drop")) %>%
    dplyr::select(
      id, 
      .pred_survival, 
      EVENT_STAT, 
      EVENT_MON, 
      .eval_time,
      BRAIN_REL
    ) 
  
  
  # 6.3.3 Create bias score for defined time
  
  outliers_bias <- outlier_summary %>%
    mutate(
      bias_score = (((1 - EVENT_STAT) - .pred_survival) ^ 2) * ((EVENT_MON - .eval_time) * ( 1 - 2 * (EVENT_STAT)))
    ) %>%
    arrange(desc(bias_score))
  
  # 6.4 Identify the highest bias patients
  
  # 6.4.1 Obtain mean and sd and then filter baed on patients higher than determined SD
  
  extreme_outliers <- 
    outliers_bias %>% 
    mutate(mean_bias = mean(bias_score),
           sd_bias = sd(bias_score)) %>% 
    filter(bias_score > (mean_bias + 2 * sd_bias))
  
  # 6.4.2 Obtain their IDs
  
  top_bias_ids <- extreme_outliers$id
  
  print(length(top_bias_ids))
  
  # 6.4.3 Identify different characteristics of patients identified as top bias
  
  print(outliers_bias %>%
          filter(id %in% top_bias_ids) %>% 
          group_by(BRAIN_REL, EVENT_STAT) %>%
          summarise(
            count = n(),
            avg_pred_event = mean(.pred_survival),
            avg_event_time = mean(EVENT_MON)
          ) %>%
          arrange(desc(count))
  )
  
  # 6.5.1 Add a column identifying patients as top bias or not
  
  outliers <- 
    outliers %>% 
    mutate(quadrant = case_when(
      id %in% top_bias_ids & EVENT_STAT == 0 ~ 2,
      id %in% top_bias_ids & EVENT_STAT == 1 ~ 1,
      TRUE ~ 0
    ))
  
  print(outliers %>% 
          group_by(EVENT_STAT) %>% 
          dplyr::count(quadrant))
  
  # 6.5.2 Plot
  
  theme_embedded <- theme_classic(base_size = 25) + 
    theme(
      legend.position = c(0.95, 0.3), # Adjust coordinates (x, y) from 0 to 1
      legend.background = element_rect(fill = alpha("white", 0.5))
    )
  
  # 6.5.2.1 Plot colored by EVENT_STAT
  
  p1 <- ggplot(outliers, aes(x = EVENT_MON, y = .pred_survival, color = factor(EVENT_STAT), shape = factor(EVENT_STAT))) +
    geom_point(size = 2, alpha = 0.7) + # Increased size and opacity
    stat_ellipse(type = "t", level = 0.95) + # Adds 95% confidence ellipse
    geom_vline(xintercept = eval_time, linetype = "dashed", color = "red") +
    scale_color_viridis_d() + 
    labs(
      title = "Event Status Distribution",
      x = paste0("Actual Event Time (Months)", eval_time),
      y = "Predicted Survival",
      color = "Event Stat",
      shape = "Event Stat"
    ) +
    theme_embedded
  
  # 6.5.2.2 Plot colored by quadrant
  
  p2 <- ggplot(outliers, aes(x = EVENT_MON, y = .pred_survival, color = factor(quadrant), shape = factor(EVENT_STAT))) +
    geom_point(size = 2, alpha = 0.7) +
    stat_ellipse(aes(group = quadrant), type = "t", level = 0.95) + 
    geom_vline(xintercept = eval_time, linetype = "dashed", color = "red") +
    scale_color_viridis_d() + 
    labs(
      title = "Quadrant Analysis",
      x = paste0("Actual Event Time (Months)", eval_time),
      y = "Predicted Survival",
      color = "Quadrant",
      shape = "Event Stat"
    ) +
    theme_embedded
  
  # 6.5.3 Combine and stack
  
  print( p1 + p2)
  
  list[[i]] <- top_bias_ids
  
}



# 7.- Other scores --------------------------------------------------------

# 7.1 Brier score

eval_results.gse2034 <- final_fit %>%
  augment(new_data = proof_genes_pt_gse2034, eval_time = c(36, 60, 120)) 

performance.gse2034 <- eval_results.gse2034 %>%
  brier_survival(truth = surv_obj, .pred)

print(performance.gse2034)

# 7.2 Martingale and Schofeild residuals

cox.zph(cox_model.gse2034)



ggcoxzph(cox.zph(cox_model.gse2034))



ggcoxdiagnostics(cox_model.gse2034, type = "martingale",
                 linear.predictions = FALSE, ggtheme = theme_bw())


