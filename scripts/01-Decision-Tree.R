# 01-Decision-Tree.R
# Purpose: End-to-end decision tree analysis of gene expression data
# Author: Seb Rauschert
# Date: 11/11/2024

#--- Load required libraries
library(tidyverse)
library(tidymodels)
library(rpart.plot)
library(patchwork)
library(viridis)
library(gt)

#--- Create directories
dirs <- c("data", "models", "results", "evaluation")
walk(dirs, ~if(!dir.exists(.x)) dir.create(.x))

#--- Set random seed for reproducibility
set.seed(123)

#--- Function to generate expression data with realistic overlap
generate_expression_data <- function(n_samples = 200, n_features = 10) {
  # Create overlapping expression patterns for each subtype
  subtype1_pattern <- function() {
    # Subtype 1 characteristics with more variation
    gene_A <- rnorm(n_samples/2, mean = 1.2, sd = 0.8)   
    gene_B <- rnorm(n_samples/2, mean = -0.9, sd = 0.7)  
    gene_C <- rnorm(n_samples/2, mean = 1.5, sd = 0.9)   
    gene_D <- rnorm(n_samples/2, mean = -1.1, sd = 0.8)  
    gene_E <- rnorm(n_samples/2, mean = 0.8, sd = 0.7)   
    
    noise_genes <- matrix(rnorm(n_samples/2 * (n_features-5), mean = 0, sd = 1.2), 
                          ncol = n_features-5)
    
    noise_matrix <- matrix(rnorm(n_samples/2 * 5, mean = 0, sd = 0.4), ncol = 5)
    genes_with_noise <- cbind(gene_A, gene_B, gene_C, gene_D, gene_E) + noise_matrix
    
    cbind(genes_with_noise, noise_genes)
  }
  
  subtype2_pattern <- function() {
    gene_A <- rnorm(n_samples/2, mean = -0.8, sd = 0.8)  
    gene_B <- rnorm(n_samples/2, mean = 0.7, sd = 0.7)   
    gene_C <- rnorm(n_samples/2, mean = -1.0, sd = 0.9)  
    gene_D <- rnorm(n_samples/2, mean = 0.9, sd = 0.8)   
    gene_E <- rnorm(n_samples/2, mean = -0.6, sd = 0.7)  
    
    noise_genes <- matrix(rnorm(n_samples/2 * (n_features-5), mean = 0, sd = 1.2), 
                          ncol = n_features-5)
    
    noise_matrix <- matrix(rnorm(n_samples/2 * 5, mean = 0, sd = 0.4), ncol = 5)
    genes_with_noise <- cbind(gene_A, gene_B, gene_C, gene_D, gene_E) + noise_matrix
    
    cbind(genes_with_noise, noise_genes)
  }
  
  # Add random fluctuation to create some ambiguous cases
  fluctuation <- function(data, probability = 0.1) {
    n_rows <- nrow(data)
    n_cols <- ncol(data)
    rows_to_change <- sample(1:n_rows, size = floor(n_rows * probability))
    for(row in rows_to_change) {
      data[row,] <- data[row,] + rnorm(n_cols, mean = 0, sd = 1.5)
    }
    return(data)
  }
  
  # Combine patterns
  expression_matrix <- rbind(subtype1_pattern(), subtype2_pattern())
  expression_matrix <- fluctuation(expression_matrix)
  
  # Create data frame
  gene_names <- c(paste0("gene_", LETTERS[1:5]), 
                  paste0("noise_gene_", 1:(n_features-5)))
  expression_df <- as.data.frame(expression_matrix)
  colnames(expression_df) <- gene_names
  expression_df$subtype <- factor(rep(c("Subtype1", "Subtype2"), each = n_samples/2))
  
  return(expression_df)
}

#--- Generate and split data
gene_expression <- generate_expression_data(n_samples = 200, n_features = 10)
gene_split <- initial_split(gene_expression, prop = 0.75, strata = subtype)
gene_train <- training(gene_split)
gene_test <- testing(gene_split)

#--- Create model specification and tune
tree_spec <- decision_tree(
  cost_complexity = tune(),
  min_n = tune(),
  tree_depth = tune()
) %>%
  set_engine("rpart") %>%
  set_mode("classification")

#--- Create cross-validation folds
cv_folds <- vfold_cv(gene_train, v = 5, strata = subtype)

# Define parameter grid
param_grid <- grid_regular(
  cost_complexity(range = c(-4, -1)),
  min_n(range = c(2, 10)),
  tree_depth(range = c(3, 10)),
  levels = 8
)

#--- Tune model
tune_results <- tune_grid(
  tree_spec,
  subtype ~ .,
  resamples = cv_folds,
  grid = param_grid,
  metrics = metric_set(accuracy, roc_auc, precision, recall)
)

#--- Select best parameters and train final model
best_params <- select_best(tune_results, metric = "roc_auc")
final_tree <- finalize_model(tree_spec, best_params) %>%
  fit(subtype ~ ., data = gene_train)

#--- Function to evaluate model performance
evaluate_performance <- function(model, data, dataset_name) {
  predictions <- predict(model, data, type = "prob") %>%
    bind_cols(predict(model, data)) %>%
    bind_cols(data %>% select(subtype)) %>%
    rename(
      prob_subtype1 = .pred_Subtype1,
      prob_subtype2 = .pred_Subtype2,
      predicted_class = .pred_class,
      actual_class = subtype
    )
  
  metrics <- predictions %>%
    metrics(truth = actual_class, 
            estimate = predicted_class, 
            prob_subtype1)
  
  conf_mat_plot <- conf_mat(predictions, 
                            truth = actual_class,
                            estimate = predicted_class) %>%
    autoplot(type = "heatmap") +
    scale_fill_viridis() +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(size = 16, face = "bold")) +
    labs(title = paste(dataset_name, "Confusion Matrix"))
  
  roc_plot <- predictions %>%
    roc_curve(truth = actual_class, prob_subtype1) %>%
    autoplot() +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(size = 16, face = "bold")) +
    labs(title = paste(dataset_name, "ROC Curve"))
  
  list(
    predictions = predictions,
    metrics = metrics,
    conf_mat_plot = conf_mat_plot,
    roc_plot = roc_plot
  )
}

#--- Create visualizations
#--- Decision tree
png("evaluation/decision_tree.png", width = 1200, height = 800, res = 120)
rpart.plot(tree_fit$fit,
           type = 5,                    # Include class probabilities
           extra = 6,                   # Show class percentages
           box.palette = "Blues",       # Blue gradient
           branch.lty = 1,              # Solid lines
           branch.lwd = 1.5,            # Medium thickness
           split.cex = 1.2,             # Split text size
           nn.font = 4,                 # Node text font (modern)
           leaf.round = 1,              # Rounded leaf nodes
           roundint = FALSE,            # Avoid rounding
           digits = 2,                  # Decimal places
           under = TRUE,                # Put text under boxes
           main = "Gene Expression Analysis")
dev.off()

#--- Evaluate performance
train_results <- evaluate_performance(final_tree, gene_train, "Training")
test_results <- evaluate_performance(final_tree, gene_test, "Test")

#--- Save performance plots
ggsave("evaluation/train_confusion_matrix.png", train_results$conf_mat_plot)
ggsave("evaluation/train_roc_curve.png", train_results$roc_plot)
ggsave("evaluation/test_confusion_matrix.png", test_results$conf_mat_plot)
ggsave("evaluation/test_roc_curve.png", test_results$roc_plot)

#--- Save metrics as CSV
bind_rows(
  train_results$metrics %>% mutate(dataset = "Training"),
  test_results$metrics %>% mutate(dataset = "Test")
) %>%
  write_csv("evaluation/model_metrics.csv")

# Save hyperparameters
write_csv(best_params, "evaluation/best_parameters.csv")