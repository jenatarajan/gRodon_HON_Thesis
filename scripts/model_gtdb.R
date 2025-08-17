
library(dplyr)
library(MASS)


#the merged data frame will contain columns: Genome, d, CUBHE, ConsistencyHE, CPB
pred_df <- readRDS("merged_codon_stats_growth.rds")
cat("Loaded data with", nrow(pred_df), "rows.\n")
cat("Columns: ", paste(names(pred_df), collapse = ", "), "\n")


#define Box–Cox Transformation Helper Function
boxcoxTransform <- function(x, lambda, back_transform = FALSE){
  if(back_transform){
    # Inverse transformation: (x * lambda + 1)^(1/lambda)
    return(((x * lambda) + 1)^(1 / lambda))
  } else {
    # Forward transformation: ((x^lambda) - 1) / lambda
    return(((x^lambda) - 1) / lambda)
  }
}

#split Data into 90% Training and 10% Testing
set.seed(123)
n <- nrow(pred_df)
train_idx <- sample(seq_len(n), size = floor(0.9 * n))
test_idx <- setdiff(seq_len(n), train_idx)
cat("Training rows:", length(train_idx), " Testing rows:", length(test_idx), "\n")

train_data <- pred_df[train_idx, ]
test_data <- pred_df[test_idx, ]

#fit Model on Training Data Using Box–Cox Transform
# Compute the optimal lambda on the training data using the model:
# d ~ CUBHE + ConsistencyHE + CPB
bc_obj <- boxcox(d ~ CUBHE + ConsistencyHE + CPB, data = train_data)
lambda <- bc_obj$x[which.max(bc_obj$y)]
cat("Optimal lambda for training data:", lambda, "\n")

#transform the response variable in training data
train_data <- train_data %>%
  mutate(y_trans = boxcoxTransform(d, lambda, back_transform = FALSE))

#fit the linear model
model <- lm(y_trans ~ CUBHE + ConsistencyHE + CPB, data = train_data)

#Predict on the Test Data and Back-Transform Predictions
preds_trans <- predict(model, newdata = test_data, interval = "confidence")
preds <- boxcoxTransform(preds_trans[,"fit"], lambda, back_transform = TRUE)

#Calculate Performance Metrics (MSE) and Log-Transform MSE
mse <- mean((preds - test_data$d)^2, na.rm = TRUE)
cat("Test MSE (raw):", mse, "\n")
log_mse <- log(mse)
cat("Test MSE (log transformed):", log_mse, "\n")

#Save Final Predictions Data Frame
final_predictions <- data.frame(
  Genome = test_data$Genome,
  Actual_d = test_data$d,
  Predicted_d = preds,
  stringsAsFactors = FALSE
)

saveRDS(final_predictions, file = "final_predictions_90_10.rds")
write.csv(final_predictions, file = "final_predictions_90_10.csv", row.names = FALSE)
cat("Final predictions saved as 'final_predictions_90_10.rds' and 'final_predictions_90_10.csv'.\n")

