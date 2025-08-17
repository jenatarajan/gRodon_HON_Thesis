
library(dplyr)
library(MASS)
library(ggplot2)
library(tidyr)

#load Your Data
load("/gpfs/projects/WeissmanGroup/gtdb220/gtdb220_genome_predictions_full.rda")
cat("Loaded full-length data with", nrow(pred_df), "rows.\n")
cat("Columns: ", paste(names(pred_df), collapse = ", "), "\n")

#the merged data frame should contain Genome, d, CUBHE, ConsistencyHE, CPB
pred_df_truncated <- readRDS("merged_codon_stats_growth.rds")
cat("Loaded 180bp data with", nrow(pred_df_truncated), "rows.\n")
cat("Columns: ", paste(names(pred_df_truncated), collapse = ", "), "\n")


#find common genomes between both data frames
common_genomes <- intersect(pred_df$Genome, pred_df_truncated$Genome)

#filter both data frames to include only rows with genomes in common
pred_df <- pred_df %>% filter(Genome %in% common_genomes)
pred_df_truncated <- pred_df_truncated %>% filter(Genome %in% common_genomes)

cat("Filtered full-length data with", nrow(pred_df), "rows.\n")
cat("Filtered 180bp data with", nrow(pred_df_truncated), "rows.\n")

#split Data into 90% Training and 10% Testing 
set.seed(123)
n <- nrow(pred_df_truncated)  # Same number of rows for both data frames after filtering
train_idx <- sample(seq_len(n), size = floor(0.9 * n))
test_idx <- setdiff(seq_len(n), train_idx)
cat("Training rows:", length(train_idx), " Testing rows:", length(test_idx), "\n")

train_data <- pred_df_truncated[train_idx, ]
test_data <- pred_df_truncated[test_idx, ]
train_data_full <- pred_df[train_idx, ]
test_data_full <- pred_df[test_idx, ]

#Define Box–Cox Transformation
boxcoxTransform <- function(x, lambda, back_transform = FALSE){
  if(back_transform){
    # Inverse transformation: (x * lambda + 1)^(1/lambda)
    return(((x * lambda) + 1)^(1 / lambda))
  } else {
    # Forward transformation: ((x^lambda) - 1) / lambda
    return(((x^lambda) - 1) / lambda)
  }
}

#fit the 180 Model on the Training Data Using Box–Cox Transform
# Compute the optimal lambda: d ~ CUBHE + ConsistencyHE + CPB
bc_obj_180 <- boxcox(d ~ CUBHE + ConsistencyHE + CPB, data = train_data)
lambda_180 <- bc_obj_180$x[which.max(bc_obj_180$y)]
cat("Optimal lambda for 180bp model training data:", lambda_180, "\n")

#transform the response variable in training data.
train_data <- train_data %>%
  mutate(y_trans_180 = boxcoxTransform(d, lambda_180, back_transform = FALSE))

#fit the linear model on the transformed response for the 180bp model
model_180 <- lm(y_trans_180 ~ CUBHE + ConsistencyHE + CPB, data = train_data)

#fit full model on training data
bc_obj_full <- boxcox(d ~ CUBHE + ConsistencyHE + CPB, data = train_data_full)
lambda_full <- bc_obj_full$x[which.max(bc_obj_full$y)]
cat("Optimal lambda for full-length model training data:", lambda_full, "\n")

#transform the response variable in training
train_data_full <- train_data_full %>%
  mutate(y_trans_full = boxcoxTransform(d, lambda_full, back_transform = FALSE))

#fit the linear model on the transformed response for the full-length model
model_full <- lm(y_trans_full ~ CUBHE + ConsistencyHE + CPB, data = train_data_full)

#predict on the Test Data and Back-Transform Predictions for Both Models
#180bp Model Predictions (on truncated test data)
preds_trans_180 <- predict(model_180, newdata = test_data, interval = "confidence")
preds_180 <- boxcoxTransform(preds_trans_180[,"fit"], lambda_180, back_transform = TRUE)

#Full-Length Model Predictions (on truncated test data)
preds_trans_full <- predict(model_full, newdata = test_data, interval = "confidence")
preds_full <- boxcoxTransform(preds_trans_full[,"fit"], lambda_full, back_transform = TRUE)

#create a New Data Frame with Actual and Predicted Values for Both Models
final_predictions <- test_data %>%
  mutate(
    Actual_d = d,  #this is the Actual d (from merged codon stats)
    Predicted_d_180 = preds_180,  #predictions from the 180bp model
    Predicted_d_full = preds_full  #predictions from the full-length model
  )


#calculate errors (actual - predicted) for both models
final_predictions <- final_predictions %>%
  mutate(
    Error_180 = Actual_d - Predicted_d_180,
    Error_full = Actual_d - Predicted_d_full
  )
cat("Columns in final_predictions: ", paste(names(final_predictions), collapse = ", "), "\n")


#remove Major Outliers Using IQR (Interquartile Range)
iqr_180 <- IQR(final_predictions$Predicted_d_180 - final_predictions$Actual_d, na.rm = TRUE)
iqr_full <- IQR(final_predictions$Predicted_d_full - final_predictions$Actual_d, na.rm = TRUE)

#set the upper and lower bounds for the errors
lower_bound_180 <- quantile(final_predictions$Predicted_d_180 - final_predictions$Actual_d, 0.25) - 1.5 * iqr_180
upper_bound_180 <- quantile(final_predictions$Predicted_d_180 - final_predictions$Actual_d, 0.75) + 1.5 * iqr_180

lower_bound_full <- quantile(final_predictions$Predicted_d_full - final_predictions$Actual_d, 0.25) - 1.5 * iqr_full
upper_bound_full <- quantile(final_predictions$Predicted_d_full - final_predictions$Actual_d, 0.75) + 1.5 * iqr_full

#filter out outliers
error_data_filtered <- final_predictions %>%
  filter((Predicted_d_180 - Actual_d >= lower_bound_180) & (Predicted_d_180 - Actual_d <= upper_bound_180)) %>%
  filter((Predicted_d_full - Actual_d >= lower_bound_full) & (Predicted_d_full - Actual_d <= upper_bound_full))

cat("Columns in final_predictions: ", paste(names(error_data_filtered), collapse = ", "), "\n")

#create a new data frame with errors
error_data <- error_data_filtered %>%
 dplyr::select(Genome, Predicted_d_180, Predicted_d_full, Actual_d) %>%
  mutate(
    Error_180 = Actual_d - Predicted_d_180,
    Error_full = Actual_d - Predicted_d_full
  ) %>%
  pivot_longer(cols = starts_with("Error"), names_to = "Model", values_to = "Error")

#create box plot for errors
p <- ggplot(error_data, aes(x = Model, y = Error, color = Model)) +
  geom_boxplot(outlier.size = 2, outlier.colour = "red", fill = "lightgray") +  # Boxplot for errors
  labs(title = "Error Distribution for Predicted Doubling Times (After Removing Outliers)",
       x = "Model", y = "Error (Actual - Predicted)") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(1, 2, 1, 2, "cm")) +
  scale_y_continuous(expand = c(0, 0.5))  

#save
pdf("Performance_Summary_Error_BoxPlot_both_models_filtered.pdf", width = 6, height = 5)
print(p)
dev.off()

cat("Boxplot after removing outliers saved as 'Performance_Summary_Error_BoxPlot_both_models_filtered.pdf' and 'Performance_Summary_Error_BoxPlot_both_models_filtered.png'.\n")


final_df <- data.frame(
  Genome = test_data$Genome,
  Actual_d = test_data$d,
  Predicted_d_180 = preds_180,
  Predicted_d_full = preds_full,
  Percent_Error_180 = abs(preds_180 - test_data$d) / abs(test_data$d) * 100,
  Percent_Error_Full = abs(preds_full - test_data$d) / abs(test_data$d) * 100
)

#save data frame as csv
write.csv(final_df, "final_predictions_with_percent_error.csv", row.names = FALSE)


plot_df <- final_df %>%
  dplyr::select(Genome, Percent_Error_180, Percent_Error_Full) %>%
  pivot_longer(cols = starts_with("Percent_Error"), 
               names_to = "Model", 
               values_to = "Percent_Error") %>%
  mutate(Model = recode(Model,
                        "Percent_Error_180" = "180bp Model",
                        "Percent_Error_Full" = "Full-Length Model"))

#save as pdf

pdf("percent_error_boxplot.pdf", width = 6, height = 5)

ggplot(plot_df, aes(x = Model, y = Percent_Error, fill = Model)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  coord_cartesian(ylim = c(0, 100)) + 
  labs(title = "Percent Error of Predictions",
       y = "Percent Error (%)",
       x = "") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

dev.off()


#calculate paired percent errors
percent_error_180 <- final_df$Percent_Error_180
percent_error_full <- final_df$Percent_Error_Full


test_result <- wilcox.test(percent_error_180, percent_error_full, paired = TRUE)
cat("Wilcoxon signed-rank test result:\n")
print(test_result)

