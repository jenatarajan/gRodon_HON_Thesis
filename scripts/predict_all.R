library(dplyr)
library(ggplot2)

#Define Gene Lengths, Directory, and Pre-Computed Lambdas
gene_lengths <- seq(30, 510, by = 30)
stats_dir <- "/gpfs/projects/WeissmanGroup/jesse_test/extdata"

#pre-computed lambdas for each gene length.
lambdas <- c(
  "30"  = -0.1010101,
  "60"  = -0.1010101,
  "90"  = -0.1010101,
  "120" = -0.1010101,
  "150" = -0.1010101,
  "180" = -0.1010101,
  "210" = -0.1414141,
  "240" = -0.1414141,
  "270" = -0.1414141,
  "300" = -0.1414141,
  "330" = -0.1414141,
  "360" = -0.1414141,
  "390" = -0.1414141,
  "420" = -0.1414141,
  "450" = -0.1414141,
  "480" = -0.1414141,
  "510" = -0.1414141
)

#Define Box-Cox Transformation Helper
boxcoxTransform <- function(x, lambda, back_transform = FALSE){
  if(back_transform){
    # Inverse transformation: (x*lambda + 1)^(1/lambda)
    return(((x * lambda) + 1)^(1 / lambda))
  } else {
    # Forward transformation: ((x^lambda) - 1) / lambda
    return(((x^lambda) - 1) / lambda)
  }
}

#generate 90%/10% Split Based on One Data Frame
#use the data from gene length 30 to define the split
split_file <- file.path(stats_dir, "StatData_Cleaned_30.rda")
if(!file.exists(split_file)) stop("File not found: ", split_file)
e_split <- new.env()
load(split_file, envir = e_split)
if(!exists("my_stat_data", envir = e_split)) stop("'my_stat_data' not found in ", split_file)
common_data <- e_split$my_stat_data
n <- nrow(common_data)
set.seed(123)
train_idx <- sample(seq_len(n), size = floor(0.9 * n))
test_idx <- setdiff(seq_len(n), train_idx)
cat("Common split generated: ", length(train_idx), " training rows and ", length(test_idx), " testing rows.\n")

#main loop: for each training length, tit model on 90% and predict on 10% of each test length
all_predictions <- data.frame()

for (train_len in gene_lengths) {
  cat("Processing training length:", train_len, "\n")
  
  #load training data for current length.
  train_file <- file.path(stats_dir, paste0("StatData_Cleaned_", train_len, ".rda"))
  if(!file.exists(train_file)) {
    cat("Train file not found for length", train_len, "- skipping.\n")
    next
  }
  e_train <- new.env()
  load(train_file, envir = e_train)
  if(!exists("my_stat_data", envir = e_train)) {
    cat("No 'my_stat_data' found in train file for length", train_len, "- skipping.\n")
    next
  }
  train_data_full <- e_train$my_stat_data
  if(is.null(train_data_full) || nrow(train_data_full) == 0) {
    cat("Train data empty for length", train_len, "- skipping.\n")
    next
  }
  
  #subset using the common split indices.
  train_subset <- train_data_full[train_idx, ]
  
  #retrieve the pre-computed lambda for this training length.
  current_lambda <- lambdas[as.character(train_len)]
  cat("Training length", train_len, ": using lambda =", current_lambda, "\n")
  
  #transform the observed doubling times in the training set
  train_subset <- train_subset %>%
    mutate(y_trans = boxcoxTransform(d, current_lambda, back_transform = FALSE))
  
  #fit a linear model on the training subset
  model <- lm(y_trans ~ CUBHE + ConsistencyHE + CPB, data = train_subset)
  
  #for each test length, load the test data and predict
  for (test_len in gene_lengths) {
    cat("  Predicting for test length:", test_len, "\n")
    test_file <- file.path(stats_dir, paste0("StatData_Cleaned_", test_len, ".rda"))
    if(!file.exists(test_file)) {
      cat("  Test file not found for length", test_len, "- skipping.\n")
      next
    }
    e_test <- new.env()
    load(test_file, envir = e_test)
    if(!exists("my_stat_data", envir = e_test)) {
      cat("  No 'my_stat_data' found in test file for length", test_len, "- skipping.\n")
      next
    }
    test_data_full <- e_test$my_stat_data
    if(is.null(test_data_full) || nrow(test_data_full) == 0) {
      cat("  Test data empty for length", test_len, "- skipping.\n")
      next
    }
    #subset test data using the test indices
    test_subset <- test_data_full[test_idx, ]
    
    #predict with box cox transformation
    preds_trans <- predict(model, newdata = test_subset, interval = "confidence")
    #back-transform predictions
    preds <- boxcoxTransform(preds_trans[,"fit"], current_lambda, back_transform = TRUE)
    lwr <- boxcoxTransform(preds_trans[,"lwr"], current_lambda, back_transform = TRUE)
    upr <- boxcoxTransform(preds_trans[,"upr"], current_lambda, back_transform = TRUE)
    
    #create a data frame of predictions
    #use "Species" from test data as the identifier
    tmp_df <- data.frame(
      TrainLength = train_len,
      TestLength = test_len,
      Species = test_subset$Species, 
      Predicted_d = preds,
      LowerCI = lwr,
      UpperCI = upr,
      Actual_d = test_subset$d,
      stringsAsFactors = FALSE
    )
    
    all_predictions <- rbind(all_predictions, tmp_df)
  }
}

#save the combined predictions
saveRDS(all_predictions, file = "all_predictions.rds")
cat("Done! Saved 'all_predictions.rds' with", nrow(all_predictions), "rows.\n")

