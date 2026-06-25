library(tidymodels)
library(censored) 
library(doParallel)
library(parallel)

# In this script we do refolds based on Monte Carlo cross validation so as to select the hyperparameters to be used in the regression

# 1.- Initial Tidymoedls workflow -----------------------------------------

# 1.1 Recipe using proof_genes_pt, so all the patients of METABRIC ER+

lr_rec <- recipe(surv_obj ~ ., data = proof_genes_pt) %>% 
  update_role(EVENT_MON, EVENT_STAT, new_role = "non_predictor") %>%
  step_dummy(all_nominal_predictors()) %>% 
  step_zv(all_predictors()) %>% 
  step_nzv(all_predictors()) %>% 
  step_normalize(all_predictors()) %>% 
  step_spatialsign()

# 1.2 Model with tune parameters

lr_mod <- proportional_hazards(
  penalty = tune(), 
  mixture = tune()
) %>%
  set_engine("glmnet")

# 1.3 Workflow 

lr_wf <- workflow() %>%
  add_model(lr_mod) %>%
  add_recipe(lr_rec)


# 2.- Define folds and tune grid ------------------------------------------


# 2.1 Define 20 resamples with Monte Carlo validation (no replacement)

set.seed(123)
resamples_20 <- mc_cv(
  proof_genes_pt, 
  prop = 0.8, 
  times = 20, 
  strata = EVENT_STAT
)

# 2.2 Grid Definition

grid <- grid_regular(
  penalty(range = c( - 4, 1)),   
  mixture(range = c(0, 1)),
  levels = 10
)


# 2.3.1 Set metrics to evaluate with tune grid levels

survival_metrics <- metric_set(
  concordance_survival,
  roc_auc_survival
)

# 2.3.2 Run grid tuning 

res_nested <- tune_grid(
  lr_wf,
  resamples = resamples_20,
  grid = grid,
  metrics = survival_metrics,
  control = control_grid(save_pred = TRUE),
  eval_time = c(36, 60, 120) 
  
)


# 3.- Tuning results ------------------------------------------------------

# 3.1 Tuning results

final_metrics <- collect_metrics(res_nested)

# 3.2 Select best parameter for c-score

best_params <- select_best(res_nested, metric = "concordance_survival")

# 3.3  Show best parameters for AUC at 3 time points

for (i in c(36, 60, 120)) {
  best_params_auc <- show_best(res_nested, metric = "roc_auc_survival", n = 20, eval_time = i)
  print(best_params_auc)
}

