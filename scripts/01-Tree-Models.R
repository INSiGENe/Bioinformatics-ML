# Load required libraries
library(tidyverse)
library(tidymodels)
library(TCGAbiolinks)
library(DESeq2)
library(xgboost)
library(vip)

#' This analysis uses breast cancer (BRCA) data from The Cancer Genome Atlas (TCGA). 
#' The data consists of RNA-seq gene expression measurements, specifically using STAR-aligned read 
#' counts, from both primary tumor samples and normal tissue samples. The preprocessing pipeline 
#' creates a balanced dataset of 30 samples (15 tumor, 15 normal) and focuses on the 100 most 
#' variable genes to reduce dimensionality while retaining informative features. The expression 
#' data is variance-stabilized using DESeq2's VST transformation to handle RNA-seq count data 
#' appropriately. This creates a normalized dataset suitable for machine learning where each row 
#' represents a sample (either tumor or normal tissue) and each column represents the expression 
#' level of one of the top 100 most variable genes. This preprocessing approach helps ensure that 
#' the machine learning models are working with standardized, biologically relevant features while
#'maintaining a balanced classification problem between tumor and normal samples.

# Load and prepare TCGA data ----------------------------------------------------
get_tcga_data <- function(n_samples = 30) {
  # Query TCGA data
  query <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = c("Primary Tumor", "Solid Tissue Normal")
  )

  GDCdownload(query)
  data <- GDCprepare(query)

  # Extract count matrix and sample info
  counts <- assay(data)
  colData <- colData(data)

  # Create balanced subset
  n_per_class <- floor(n_samples/2)
  normal_idx <- which(colData$sample_type == "Solid Tissue Normal")[1:n_per_class]
  tumor_idx <- which(colData$sample_type == "Primary Tumor")[1:n_per_class]
  keep_idx <- c(normal_idx, tumor_idx)

  # Create DESeqDataSet
  subset_counts <- counts[, keep_idx]
  subset_colData <- colData[keep_idx,]
  dds <- DESeqDataSetFromMatrix(
    countData = subset_counts,
    colData = subset_colData,
    design = ~ sample_type
  )

  # Transform data
  vsd <- vst(dds, blind = FALSE)
  expr_matrix <- assay(vsd)

  # Calculate variance for each gene and select top 100 most variable
  gene_vars <- apply(expr_matrix, 1, var)
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:100]

  # Subset expression matrix to top variable genes
  expr_matrix_filtered <- expr_matrix[top_genes, ]

  # Prepare classification data
  sample_data <- t(expr_matrix_filtered) %>%
    as.data.frame() %>%
    mutate(treatment = factor(subset_colData$sample_type))

  return(sample_data)
}

# Get data and prepare splits
set.seed(123)
sample_data <- get_tcga_data(30)

qs::qsave(sample_data, "data/sample_data.qs")

sample_data <- qs::qread("data/sample_data.qs")

# Initial split into train+validation and test sets
splits <- initial_split(sample_data, prop = 0.75, strata = treatment)
train_data <- training(splits)
test_data <- testing(splits)

# Create validation set from training data
set.seed(234)
val_set <- validation_split(train_data, 
                            prop = 0.80,
                            strata = treatment)

# Create preprocessing recipe ---------------------------------------------
rna_recipe <- recipe(treatment ~ ., data = train_data) %>%
  step_normalize(all_predictors()) %>%
  step_log(all_predictors(), offset = 1e-6)

# Model Specifications -------------------------------------------------
# 1. Decision Tree
tree_spec <- decision_tree(
  cost_complexity = tune()
) %>%
  set_engine("rpart") %>%
  set_mode("classification")

tree_workflow <- workflow() %>%
  add_recipe(rna_recipe) %>%
  add_model(tree_spec)

tree_grid <- tibble(cost_complexity = c(0.01, 0.001, 0.0001))

# 2. Random Forest
rf_spec <- rand_forest(
  mtry = tune(),
  trees = 500
) %>%
  set_engine("ranger", importance = "impurity") %>%
  set_mode("classification")

rf_workflow <- workflow() %>%
  add_recipe(rna_recipe) %>%
  add_model(rf_spec)

rf_grid <- tibble(mtry = c(10, 20, 30))

# 3. XGBoost
xgb_spec <- boost_tree(
  trees = 100,
  tree_depth = tune()
) %>%
  set_engine("xgboost") %>%
  set_mode("classification")

xgb_workflow <- workflow() %>%
  add_recipe(rna_recipe) %>%
  add_model(xgb_spec)

xgb_grid <- tibble(tree_depth = c(3, 4, 6))

# Model Evaluation Function --------------------------------------------
fit_and_evaluate <- function(workflow, grid, model_name) {
  tune_results <- workflow %>%
    tune_grid(
      val_set,
      grid = grid,
      control = control_grid(save_pred = TRUE),
      metrics = metric_set(roc_auc, accuracy)
    )
  
  best_params <- select_best(tune_results, metric = "roc_auc")
  
  levels_treatment <- levels(train_data$treatment)
  target_level <- levels_treatment[1]
  
  roc_data <- tune_results %>%
    collect_predictions(parameters = best_params) %>%
    roc_curve(treatment, paste0(".pred_", target_level)) %>%
    mutate(model = model_name)
  
  list(
    tune_results = tune_results,
    best_params = best_params,
    roc_data = roc_data
  )
}

# Evaluate Models ----------------------------------------------------
tree_results <- fit_and_evaluate(tree_workflow, tree_grid, "Decision Tree")
rf_results <- fit_and_evaluate(rf_workflow, rf_grid, "Random Forest")
xgb_results <- fit_and_evaluate(xgb_workflow, xgb_grid, "XGBoost")

# Visualization Function ---------------------------------------------
create_model_plots <- function(rf_results, final_fit, predictions) {
  # 1. Variable Importance Plot
  importance_plot <- final_fit %>%
    extract_fit_parsnip() %>%
    vip::vi() %>%
    mutate(Variable = fct_reorder(Variable, Importance)) %>%
    slice_max(Importance, n = 10) %>%
    ggplot(aes(x = Variable, y = Importance)) +
    geom_col(aes(fill = Importance)) +
    geom_text(aes(label = sprintf("%.3f", Importance)), hjust = -0.1) +
    scale_fill_viridis_c() +
    coord_flip() +
    theme_minimal() +
    labs(title = "Top 10 Most Important Genes")
  
  # 2. ROC Curve
  roc_curve_plot <- predictions %>%
    roc_curve(truth = treatment, 
              ".pred_Primary Tumor") %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity)) +
    geom_path(size = 1.2) +
    geom_abline(lty = 3) +
    coord_equal() +
    theme_minimal() +
    labs(title = "ROC Curve on Test Set")
  
  # 3. Confusion Matrix
  conf_mat_plot <- predictions %>%
    conf_mat(truth = treatment, estimate = .pred_class) %>%
    autoplot(type = "heatmap") +
    theme_minimal() +
    labs(title = "Confusion Matrix")
  
  # Combine plots using patchwork
  library(patchwork)
  combined_plot <- (importance_plot + roc_curve_plot) / conf_mat_plot
  
  return(combined_plot)
}

# Final Model Fitting and Visualization -------------------------------
# Function to fit final model and get predictions
fit_final_model <- function(model_spec, model_results, splits, recipe) {
  final_workflow <- workflow() %>%
    add_recipe(recipe) %>%
    add_model(model_spec)
  
  final_fit <- final_workflow %>%
    last_fit(splits)
  
  return(final_fit)
}

# Fit final Decision Tree model
best_tree <- decision_tree(
  cost_complexity = tree_results$best_params$cost_complexity
) %>%
  set_engine("rpart") %>%
  set_mode("classification")

final_tree_fit <- fit_final_model(best_tree, tree_results, splits, rna_recipe)

# Fit final Random Forest model
best_rf <- rand_forest(
  mtry = rf_results$best_params$mtry,
  trees = 500
) %>%
  set_engine("ranger", importance = "impurity") %>%
  set_mode("classification")

final_rf_fit <- fit_final_model(best_rf, rf_results, splits, rna_recipe)

# Fit final XGBoost model
best_xgb <- boost_tree(
  trees = 100,
  tree_depth = xgb_results$best_params$tree_depth
) %>%
  set_engine("xgboost") %>%
  set_mode("classification")

final_xgb_fit <- fit_final_model(best_xgb, xgb_results, splits, rna_recipe)

# Create plots for each model
create_model_plots <- function(final_fit, predictions, model_name) {
  # 1. Variable Importance Plot
  importance_plot <- final_fit %>%
    extract_fit_parsnip() %>%
    vip::vi() %>%
    mutate(Variable = fct_reorder(Variable, Importance)) %>%
    slice_max(Importance, n = 10) %>%
    ggplot(aes(x = Variable, y = Importance)) +
    geom_col(aes(fill = Importance)) +
    geom_text(aes(label = sprintf("%.3f", Importance)), hjust = -0.1) +
    scale_fill_viridis_c() +
    coord_flip() +
    theme_minimal() +
    labs(title = paste("Top 10 Most Important Genes -", model_name))
  
  # 2. ROC Curve
  roc_curve_plot <- predictions %>%
    roc_curve(truth = treatment, 
              ".pred_Primary Tumor") %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity)) +
    geom_path(size = 1.2) +
    geom_abline(lty = 3) +
    coord_equal() +
    theme_minimal() +
    labs(title = paste("ROC Curve on Test Set -", model_name))
  
  # 3. Confusion Matrix
  conf_mat_plot <- predictions %>%
    conf_mat(truth = treatment, estimate = .pred_class) %>%
    autoplot(type = "heatmap") +
    theme_minimal() +
    labs(title = paste("Confusion Matrix -", model_name))
  
  # Combine plots using patchwork
  combined_plot <- (importance_plot + roc_curve_plot) / conf_mat_plot
  
  return(combined_plot)
}

# Get predictions and create visualizations for each model
tree_predictions <- collect_predictions(final_tree_fit)
rf_predictions <- collect_predictions(final_rf_fit)
xgb_predictions <- collect_predictions(final_xgb_fit)

# Create plots for each model
tree_plots <- create_model_plots(final_tree_fit, tree_predictions, "Decision Tree")
rf_plots <- create_model_plots(final_rf_fit, rf_predictions, "Random Forest")
xgb_plots <- create_model_plots(final_xgb_fit, xgb_predictions, "XGBoost")

# Print final performance metrics for each model
print("Decision Tree Test Set Performance:")
print(collect_metrics(final_tree_fit))

print("Random Forest Test Set Performance:")
print(collect_metrics(final_rf_fit))

print("XGBoost Test Set Performance:")
print(collect_metrics(final_xgb_fit))

# Display all plots
# You can view them individually:
print(tree_plots)
print(rf_plots)
print(xgb_plots)

sample_data |> ggplot(aes(y = ENSG00000168079.17, x = treatment, fill = treatment)) + geom_boxplot(alpha=0.8) + theme_bw()