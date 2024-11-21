# Tree based Models from Start to Finish: A Beginner's Guide with R

In this tutorial, we'll walk through a practical example of applying machine learning to biological data. We'll use breast cancer gene expression data as our example dataset - not because we're trying to make medical discoveries, but because it's a great way to learn how to work with high-dimensional data in a real-world context.

## Understanding Our Example Dataset

We're using publicly available data from The Cancer Genome Atlas (TCGA) as our learning dataset. This kind of data is perfect for learning machine learning concepts because:

- It has a clear binary classification goal (distinguishing between two types of samples)
- It contains many features (genes)
- It represents real-world complexity
- It's freely available for educational purposes

For our learning exercise, we're working with 30 samples (15 of each type) and 100 genes. We've intentionally kept the dataset small and manageable - perfect for understanding the fundamental concepts without getting lost in complexity.

## Machine Learning Models: A Conceptual Overview

Before diving into the code, let's understand the three different types of models we'll be using:

### Decision Trees

Think of these as a series of yes/no questions about gene expression levels. For example, "Is Gene A's expression above X?" If yes, go one way; if no, go another way. Keep asking questions until you reach a final decision about whether the sample is tumor or normal. They're like a flowchart that makes decisions based on data.

Key advantages:

- Easy to understand and visualize
- Can handle different types of data
- Make intuitive sense

### Random Forests

This is like having a committee of decision trees, each looking at different combinations of genes. They all vote on whether a sample is tumor or normal, and we take the majority decision. It's like having multiple experts give their opinion and taking a vote.

Key advantages:

- More robust than single decision trees
- Can handle complex relationships in data
- Less likely to overfit (memorize) the training data

### XGBoost

This is a more sophisticated version of the committee approach, where each new tree focuses specifically on the mistakes made by previous trees. It's like having each new committee member specialize in the cases where previous members struggled.

Key advantages:

- Often provides better predictions than simpler models
- Automatically handles complex interactions
- Very popular in real-world applications

## Introduction to Our Tools

Before we dive in, let's understand the key libraries we're using:

```r
library(tidyverse)    # A collection of R packages for data manipulation and visualization
library(tidymodels)   # A collection of packages for modeling and machine learning
library(TCGAbiolinks) # For accessing cancer genomics data
library(DESeq2)       # For processing gene expression data
library(xgboost)      # For gradient boosting machines
library(vip)          # For variable importance plots

```

Let's break down what each library does:

- `tidyverse`: Think of this as your Swiss Army knife for data science. It includes tools for reading data (`readr`), manipulating data (`dplyr`), and creating visualizations (`ggplot2`).
- `tidymodels`: This is our machine learning framework. It makes complex modeling tasks more consistent and easier to understand.
- The others are specialized tools for working with biological data and creating visualizations.

## Data Preparation: A Deep Dive

The first step in any machine learning project is getting your data ready. Let's break down our data preparation function:

```r
get_tcga_data <- function(n_samples = 30) {
  # Query TCGA data
  query <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = c("Primary Tumor", "Solid Tissue Normal")
  )

```

This function is doing several important things:

1. We're creating a function that will get us 30 samples by default
2. We're querying a database (TCGA) for specific types of data
3. The `c("Primary Tumor", "Solid Tissue Normal")` part tells R we want both types of samples

Next, we process this data:

```r
  # Create balanced subset
  n_per_class <- floor(n_samples/2)  # Divide our total desired samples by 2
  normal_idx <- which(colData$sample_type == "Solid Tissue Normal")[1:n_per_class]
  tumor_idx <- which(colData$sample_type == "Primary Tumor")[1:n_per_class]
  keep_idx <- c(normal_idx, tumor_idx)

```

This code ensures we have an equal number of each type of sample - this is called "balancing" your dataset. It's like making sure you have the same number of heads and tails when teaching someone about coin flips.

### Data Transformation

```r
  # Transform data
  vsd <- vst(dds, blind = FALSE)
  expr_matrix <- assay(vsd)

  # Calculate variance for each gene and select top 100 most variable
  gene_vars <- apply(expr_matrix, 1, var)
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:100]

```

This section is crucial:

1. First, we transform our data to make it more suitable for analysis (think of it like converting temperatures from Fahrenheit to Celsius)
2. Then we calculate how much each gene varies across samples
3. We keep only the 100 genes that show the most variation - this is called "feature selection"

## Setting Up Our Analysis Pipeline

Now comes the exciting part - preparing our data for machine learning:

```r
# Initial split into train+validation and test sets
splits <- initial_split(sample_data, prop = 0.75, strata = treatment)
train_data <- training(splits)
test_data <- testing(splits)

```

This is like dividing students into groups:

- Training data (75%): This is what we use to teach our model
- Testing data (25%): This is what we use to give our model its final exam

We also create a validation set:

```r
# Create validation set from training data
val_set <- validation_split(train_data,
                          prop = 0.80,
                          strata = treatment)

```

Think of the validation set as practice tests before the final exam.

### Creating Our Recipe

```r
# Create preprocessing recipe
rna_recipe <- recipe(treatment ~ ., data = train_data) %>%
  step_normalize(all_predictors()) %>%
  step_log(all_predictors(), offset = 1e-6)

```

This is like writing down the steps for a cooking recipe:

1. `treatment ~ .` means "use all columns to predict the treatment column"
2. `step_normalize` scales all our numbers to be comparable (like converting different currencies to dollars)
3. `step_log` helps handle very large numbers (like using scientific notation)

## Understanding Our Models

### Decision Trees

```r
# Decision Tree Specification
tree_spec <- decision_tree(
  cost_complexity = tune()
) %>%
  set_engine("rpart") %>%
  set_mode("classification")

```

A decision tree is like playing 20 questions:

- The model asks questions about our data ("Is this gene's expression above X?")
- Based on the answer, it asks another question
- Eventually, it makes a prediction
- `cost_complexity` controls how complex our tree can get (like limiting the number of questions)

### Random Forests

```r
# Random Forest Specification
rf_spec <- rand_forest(
  mtry = tune(),
  trees = 500
) %>%
  set_engine("ranger", importance = "impurity") %>%
  set_mode("classification")

```

Random forests are like having multiple decision trees work together:

- We create 500 different trees (`trees = 500`)
- Each tree looks at a random subset of genes (`mtry = tune()`)
- They vote on the final prediction
- It's like asking 500 different experts and taking their collective opinion

### XGBoost

```r
# XGBoost Specification
xgb_spec <- boost_tree(
  trees = 100,
  tree_depth = tune()
) %>%
  set_engine("xgboost") %>%
  set_mode("classification")

```

XGBoost is like an advanced version of random forests:

- Each new tree focuses on correcting the mistakes of previous trees
- `tree_depth` controls how complex each individual tree can be
- It's like having each new expert specifically study cases where previous experts made mistakes

## Evaluating Our Models

Here's how we evaluate each model:

```r
# Model Evaluation Function
fit_and_evaluate <- function(workflow, grid, model_name) {
  tune_results <- workflow %>%
    tune_grid(
      val_set,
      grid = grid,
      control = control_grid(save_pred = TRUE),
      metrics = metric_set(roc_auc, accuracy)
    )

```

This function:

1. Takes a model workflow (our recipe + model)
2. Tries different combinations of model settings (like adjusting the temperature when cooking)
3. Measures how well each combination performs
4. Saves the results for comparison

### Creating Visualizations

```r
# Visualization Function
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
    labs(title = "Top 10 Most Important Features")

```

This creates three important visualizations:

1. Variable Importance Plot: Shows which genes were most helpful in making predictions
2. ROC Curve: Shows how well our model balances between different types of errors
3. Confusion Matrix: Shows exactly what kinds of correct and incorrect predictions we made

## Final Model and Predictions

```r
# Final Model Fitting
best_rf <- rand_forest(
  mtry = rf_results$best_params$mtry,
  trees = 500
) %>%
  set_engine("ranger", importance = "impurity") %>%
  set_mode("classification")

final_workflow <- workflow() %>%
  add_recipe(rna_recipe) %>%
  add_model(best_rf)

final_fit <- final_workflow %>%
  last_fit(splits)

```

This final section:

1. Takes our best-performing model settings
2. Creates a final model with those settings
3. Fits it to our training data
4. Makes predictions on our test data (the final exam!)

## Understanding the Results

When we run this analysis, we get several key pieces of information:

1. How accurate our model is
2. Which genes were most important for making predictions
3. Visual confirmation of our model's performance

The beauty of this approach is that while we're using gene expression data here, the same code structure could work for many other types of data:

- Customer purchase history
- Weather predictions
- Image classification
- Text analysis

## Key Takeaways for Beginners

1. Data preparation is crucial and often takes the most time
2. Always split your data into training and testing sets
3. Create a clear recipe for processing your data
4. Try multiple types of models - they each have strengths
5. Visualize your results to understand what's happening
6. Use validation to make sure your model will work on new data

Remember, this is a learning example. The same principles apply whether you're analyzing genes, predicting stock prices, or categorizing cat photos!