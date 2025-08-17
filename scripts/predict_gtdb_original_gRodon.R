library(gRodon, lib.loc="/gpfs/projects/WeissmanGroup/jesse_test")
library(Biostrings)
library(dplyr)

# Specify the text file that contains the list of file paths.
# Each line in the text file should be one file path.
file_paths <- readLines("/gpfs/projects/WeissmanGroup/gtdb220/ffn2.files")
setwd("/gpfs/projects/WeissmanGroup/gtdb220")
# Initialize a list to store results
results_list <- list()

# Loop over each file path
for(fp in file_paths) {
  cat("Processing file:", fp, "\n")
  
  # Load the genome sequences from the file using readDNAStringSet
  genes <- tryCatch({
    readDNAStringSet(fp)
  }, error = function(e) {
    cat("Error reading", fp, ":", e$message, "\n")
    return(NULL)
  })
  
  if(is.null(genes) || length(genes) == 0) {
    cat("No sequences found in", fp, "\n")
    next
  }
  
  # Identify highly expressed genes 
  highly_expressed <- grepl("ribosomal protein", names(genes), ignore.case = TRUE)
  
  # Check if there are at least 10 highly expressed genes
  if(sum(highly_expressed) < 10) {
    cat("Not enough highly expressed genes in", fp, "\n")
    next
  }
  
  # Run the gRodon prediction.
  prediction <- tryCatch({
    predictGrowth(genes, highly_expressed)
  }, error = function(e) {
    cat("Error predicting growth for", fp, ":", e$message, "\n")
    return(NULL)
  })
  
  if(is.null(prediction)) next
  
  # Extract predicted doubling time from the prediction list.
  predicted_d <- prediction$d
  
  # Create a data frame row with the genome name and predicted doubling time.
  genome_name <- basename(fp)
  result_row <- data.frame(
    Genome = genome_name,
    Predicted_d = predicted_d,
    stringsAsFactors = FALSE
  )
  
  results_list[[length(results_list) + 1]] <- result_row
}

# Combine the list of results into a single data frame
results_df <- bind_rows(results_list)

# Save the results to an RDS file and a CSV file
saveRDS(results_df, file = "gRodon_predictions.rds")
write.csv(results_df, file = "gRodon_predictions.csv", row.names = FALSE)

cat("Done! Processed", nrow(results_df), "genomes.\n")

