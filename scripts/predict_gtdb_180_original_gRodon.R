#Predict with gtdb genomes that have all genes cut to 180 base pairs, using original gRodon model
library(Biostrings)
library(gRodon, lib.loc="/gpfs/projects/WeissmanGroup/jesse_test")
library(dplyr)
library(parallel)

#get input directory and list all .ffn files
input_dir <- "/gpfs/projects/WeissmanGroup/gtdb220/prokka_out_files_180"
file_paths <- list.files(input_dir, pattern = "\\.ffn$", full.names = TRUE)
cat("Found", length(file_paths), "files in", input_dir, "\n")

#define a function to process each genome file
predictGenome <- function(fp) {
  cat("Processing file:", fp, "\n")
  
  #read the genome sequences from the file
  genes <- tryCatch({
    readDNAStringSet(fp)
  }, error = function(e) {
    cat("Error reading file:", fp, "\n", e$message, "\n")
    return(NULL)
  })
  
  #skip if no sequences found
  if (is.null(genes) || length(genes) == 0) return(NULL)
  
  #result is a DNAStringSet
  genes <- as(genes, "DNAStringSet")
  
  #identify highly expressed genes
  highly_expressed <- grepl("ribosomal protein", names(genes), ignore.case = TRUE)
  
  #ensure there are at least 10 highly expressed genes, otherwise skip this genome.
  if (sum(highly_expressed) < 10) {
    cat("Not enough highly expressed genes in:", fp, "\n")
    return(NULL)
  }
  
  # Run gRodon prediction using predictGrowth()
  growth <- tryCatch({
    predictGrowth(genes, highly_expressed, mode = "full", training_set = "madin", fragments = TRUE)
  }, error = function(e) {
    cat("Error predicting growth for:", fp, "\n", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(growth)) return(NULL)
  
  # Build a data frame with the predictions.
  # We'll use the file's basename as the genome identifier.
  result <- data.frame(
    Genome = basename(fp),
    CUBHE = growth$CUBHE,
    CUB = growth$CUB,
    dCUB = if (!is.null(growth$CUB) && !is.null(growth$CUBHE)) (growth$CUB - growth$CUBHE) / growth$CUB else NA,
    CPB = growth$CPB,
    ConsistencyHE = growth$ConsistencyHE,
    LowerCI = growth$LowerCI,
    UpperCI = growth$UpperCI,
    d = growth$d,
    nHE = growth$nHE,
    nGenes = length(genes),
    stringsAsFactors = FALSE
  )
  
  return(result)
}


#run predictions in parallel using mclapply
pred_list <- mclapply(file_paths, predictGenome, mc.cores = 28, mc.preschedule = FALSE)
pred_df_180 <- do.call("rbind", pred_list)

#save the Predictions Data Frame
output_rds <- "/gpfs/projects/WeissmanGroup/gtdb220/genome_predictions_truncated.rds"
output_csv <- "/gpfs/projects/WeissmanGroup/gtdb220/genome_predictions_truncated.csv"
saveRDS(pred_df_180, file = output_rds)
write.csv(pred_df_180, file = output_csv, row.names = FALSE)
cat("Done! Processed", nrow(pred_df), "genomes.\n")

