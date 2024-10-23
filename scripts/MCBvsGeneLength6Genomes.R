
library(coRdon)
library(Biostrings)
library(ggplot2)
library(tidyr)
library(dplyr)

# Paths to the 6 genome files
genome_paths <- c(
  "/Users/jessenatarajan/Downloads/ncbi_dataset (11)/ncbi_dataset/data/GCF_000005845.2/cds_from_genomic.fna",
  "/Users/jessenatarajan/Downloads/ncbi_dataset (11)/ncbi_dataset/data/GCF_000008865.2/cds_from_genomic.fna",
  "/Users/jessenatarajan/Downloads/ncbi_dataset (11)/ncbi_dataset/data/GCF_000010385.1/cds_from_genomic.fna",
  "/Users/jessenatarajan/Downloads/ncbi_dataset (11)/ncbi_dataset/data/GCF_000013265.1/cds_from_genomic.fna",
  "/Users/jessenatarajan/Downloads/ncbi_dataset (11)/ncbi_dataset/data/GCF_000019645.1/cds_from_genomic.fna",
  "/Users/jessenatarajan/Downloads/ncbi_dataset (11)/ncbi_dataset/data/GCF_003697165.2/cds_from_genomic.fna"
)

# Corresponding genome IDs for labeling
genome_ids <- c(
  "GCF_000005845.2",
  "GCF_000008865.2",
  "GCF_000010385.1",
  "GCF_000013265.1",
  "GCF_000019645.1",
  "GCF_003697165.2"
)

# Define the range of gene lengths
gene_lengths <- seq(3, 510, by = 3)

# Initialize a data frame to store all results
all_results_df <- data.frame()

# Iterate over each genome
for (i in seq_along(genome_paths)) {
  # Load the genome data
  genes <- readDNAStringSet(genome_paths[i])

  # Initialize a data frame for this genome
  results_df <- data.frame()

  # Iterate over gene lengths
  for (length in gene_lengths) {
    # Truncate genes to current length
    truncated_genes <- narrow(genes[width(genes) >= length], start = 1, end = length)

    if (length(truncated_genes) > 0) {
      # Create codon usage object
      codon_counts <- codonTable(truncated_genes)

      # Calculate MILC metric
      mcb_values <- MCB(codon_counts)

      # Compute average MILC
      avg_mcb <- mean(mcb_values, na.rm = TRUE)

      # Store results
      temp_df <- data.frame(
        GeneLength = length,
        MetricValue = avg_mcb,
        GenomeID = genome_ids[i]
      )
      results_df <- rbind(results_df, temp_df)
    }
  }

  # Combine with all results
  all_results_df <- rbind(all_results_df, results_df)
}

# Plot the results
plot <- ggplot(all_results_df, aes(x = GeneLength, y = MetricValue, color = GenomeID)) +
  geom_line() +
  labs(title = "MCB vs Gene Length across Genomes",
       x = "Gene Length (Base Pairs)",
       y = "Average MCB Value") +
  theme_minimal()

# Display the plot
print(plot)

# Save the plot
#ggsave("ENC_vs_GeneLength_AllGenomes.png", plot = plot, width = 10, height = 6)
