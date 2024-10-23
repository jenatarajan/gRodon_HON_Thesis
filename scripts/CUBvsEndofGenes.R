library(gRodon)
library(Biostrings)
library(ggplot2)

# Paths to the 6 genome files (replace these with actual paths to your downloaded genomes)
genome_paths <- c(
  "/Users/jessenatarajan/Downloads/ncbi_dataset (11)/ncbi_dataset/data/GCF_000005845.2/cds_from_genomic.fna",
  "/Users/jessenatarajan/Downloads/ncbi_dataset (11)/ncbi_dataset/data/GCF_000008865.2/cds_from_genomic.fna",
  "/Users/jessenatarajan/Downloads/ncbi_dataset (11)/ncbi_dataset/data/GCF_000010385.1/cds_from_genomic.fna",
  "/Users/jessenatarajan/Downloads/ncbi_dataset (11)/ncbi_dataset/data/GCF_000013265.1/cds_from_genomic.fna",
  "/Users/jessenatarajan/Downloads/ncbi_dataset (11)/ncbi_dataset/data/GCF_000019645.1/cds_from_genomic.fna",
  "/Users/jessenatarajan/Downloads/ncbi_dataset (11)/ncbi_dataset/data/GCF_003697165.2/cds_from_genomic.fna"
)

# Define the range of gene lengths (starting from 3 to 510, incrementing by 3)
gene_lengths <- seq(3, 510, by = 3)

# Initialize an empty data frame to store CUB values, gene lengths, and genome IDs
plot_data <- data.frame()

# Iterate through each genome
for (i in seq_along(genome_paths)) {
  # Load the genome data
  genes <- readDNAStringSet(genome_paths[i])

  # Identify highly expressed genes
  highly_expressed <- grepl("ribosomal protein", names(genes), ignore.case = TRUE)

  # Initialize temporary lists to store CUB and length values for this genome
  cub_values <- c()
  length_values <- c()

  # Iterate over the range of gene lengths
  for (length in gene_lengths) {
    # Truncate the genes by cutting from the front
    truncated_genes <- narrow(genes[width(genes) >= length],
                              start = width(genes[width(genes) >= length]) - length + 1,
                              end = width(genes[width(genes) >= length]))

    # Check if there are any genes of this length to process
    if (length(truncated_genes) > 0) {
      # Subset the highly expressed vector to match the truncated genes
      highly_expressed_subset <- highly_expressed[width(genes) >= length]

      # Only run gRodon prediction if there are at least 10 highly expressed genes
      if (sum(highly_expressed_subset) >= 10) {
        # Try to run the gRodon prediction, handle errors gracefully
        tryCatch({
          results <- predictGrowth(truncated_genes, highly_expressed_subset, fragments = TRUE)

          # Check if results contain valid CUB values
          if (!is.null(results$CUB) && all(!is.na(results$CUB))) {
            # Extract the CUB (MILC) from the results
            cub_values <- c(cub_values, results$CUBHE)
            length_values <- c(length_values, rep(length, length(results$CUBHE)))
          } else {
            warning(paste("Invalid CUB values for length", length, "in genome", i))
          }
        }, error = function(e) {
          warning(paste("Error processing length", length, "in genome", i, ":", e$message))
        })
      } else {
        warning(paste("Skipping length", length, "for genome", i, ": less than 10 highly expressed genes"))
      }
    }
  }

  # Check if cub_values and length_values are not empty before creating a data frame
  if (length(cub_values) > 0 && length(length_values) > 0) {
    # Combine results for this genome with the genome ID
    genome_data <- data.frame(GeneLength = length_values, CUB = cub_values, GenomeID = paste("Genome", i))

    # Add the genome-specific data to the overall plot data
    plot_data <- rbind(plot_data, genome_data)
  } else {
    warning(paste("No data for genome", i))
  }
}

# Plot CUB values against gene length for each genome, using different lines/colors
ggplot(plot_data, aes(x = GeneLength, y = CUB, color = GenomeID, group = GenomeID)) +
  geom_line(alpha = 0.6) +
  labs(title = "CUBHE (MILC) vs Gene Length for 6 Genomes",
       x = "Gene Length (Base Pairs)",
       y = "CUBHE (MILC)") +
  theme_minimal() +
  ylim(0.2, 1) +
  theme(legend.title = element_blank())


