
library(dplyr)
library(ggplot2)
.libPaths("/gpfs/projects/WeissmanGroup/jesse_test")
install.packages("viridis", repos = "https://cran.r-project.org", lib = "/gpfs/projects/WeissmanGroup/jesse_test")

library(viridis) 

#load the merged fata frame
df_merged <- readRDS("all_predictions_327.rds")

cat("Loaded merged data. Dimensions:\n")
print(dim(df_merged))
cat("Sample rows:\n")
print(head(df_merged))


#compute Performance Grid (MSE between Predicted and Actual Doubling Times)
threshold <- 1e5
df_merged_filtered <- df_merged %>%
  mutate(sq_error = (Predicted_d - Actual_d)^2) %>%
  filter(sq_error < threshold)


perf_grid <- df_merged_filtered %>%
  filter(TestLength > 150) %>%     # ignore test lengths less than 150
  group_by(TrainLength, TestLength) %>%
  summarise(
    n = n(),
    mse = mean((Predicted_d - Actual_d)^2, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nPerformance table:\n")
print(perf_grid)
  save(df_merged_filtered, file = "PerformanceData.rda")

#plot the Performance Grid
p_mse <- ggplot(perf_grid, aes(x = factor(TrainLength), y = factor(TestLength), fill = mse)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log10", option = "viridis", direction = 1, na.value = "grey80", name = "MSE") +
  labs(
    title = "Performance Grid: MSE of Predicted vs. Actual Doubling Times",
    x = "Training Length (bp)",
    y = "Testing Length (bp)"
  ) +
  theme_minimal()

#Save the Plot
pdf("Performance_Grid_MSE_log.pdf", width = 6, height = 5)
print(p_mse)
dev.off()

png("Performance_Grid_MSE.png", width = 6, height = 5, units = "in", res = 300, type = "cairo")
print(p_mse)
dev.off()

