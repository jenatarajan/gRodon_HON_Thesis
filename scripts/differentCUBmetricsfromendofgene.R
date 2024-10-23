# Install and load necessary packages
install.packages("coRdon")
library(coRdon)
library(gRodon)
library(Biostrings)
library(ggplot2)
library(tidyr)

# Load the genome data (using an example dataset from coRdon)
path_to_genome <- system.file("extdata", "GCF_000349925.2_ASM34992v2_cds_from_genomic.fna.gz", package = "gRodon")
genes <- readDNAStringSet(path_to_genome)

# Define the range of gene lengths (starting from 3 to 510, incrementing by 3)
gene_lengths <- seq(3, 510, by = 3)

# Initialize a data frame to store results
results_df <- data.frame()

# Iterate over the range of gene lengths
for (length in gene_lengths) {
  # Truncate the genes to the current length
  truncated_genes <- narrow(genes[width(genes) >= length], start = 1, end = length)

  # Check if there are any genes of this length to process
  if (length(truncated_genes) > 0) {
    # Create a codon usage object
    codon_counts <- codonTable(truncated_genes)

    # Calculate CUB metrics
    milc_values <- MILC(codon_counts)
    enc_values <- ENC(codon_counts)
    enc_prime_values <- ENCprime(codon_counts)
    scuo_values <- SCUO(codon_counts)
    b_values <- B(codon_counts)
    mcb_values <- MCB(codon_counts)

    # For each metric, compute the average value across the genes
    avg_milc <- mean(milc_values, na.rm = TRUE)
    avg_enc <- mean(enc_values, na.rm = TRUE)
    avg_enc_prime <- mean(enc_prime_values, na.rm = TRUE)
    avg_scuo <- mean(scuo_values, na.rm = TRUE)
    avg_b <- mean(b_values, na.rm = TRUE)
    avg_mcb <- mean(mcb_values, na.rm = TRUE)

    # Append the results to the data frame
    temp_df <- data.frame(
      GeneLength = length,
      MILC = avg_milc,
      ENC = avg_enc,
      ENC_Prime = avg_enc_prime,
      SCUO = avg_scuo,
      B = avg_b,
      MCB = avg_mcb
    )
    results_df <- rbind(results_df, temp_df)
  }
}

# Reshape the data frame to long format for plotting
plot_data <- gather(results_df, key = "Metric", value = "Value", -GeneLength)

# Plot the metrics against gene length
ggplot(plot_data, aes(x = GeneLength, y = Value, color = Metric)) +
  geom_line() +
  labs(title = "CUB Metrics vs Gene Length",
       x = "Gene Length (Base Pairs)",
       y = "Average CUB Metric Value") +
  theme_minimal()

