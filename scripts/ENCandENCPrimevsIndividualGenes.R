# Install and load necessary packages
if (!require("coRdon")) install.packages("coRdon", dependencies = TRUE)
if (!require("Biostrings")) BiocManager::install("Biostrings")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("hexbin")) install.packages("hexbin")

library(coRdon)
library(Biostrings)
library(ggplot2)
library(dplyr)
library(hexbin)

# Path to the genome file (update as needed)
path_to_genome <- system.file("extdata", "GCF_000349925.2_ASM34992v2_cds_from_genomic.fna.gz", package = "gRodon")
genes <- readDNAStringSet(path_to_genome)

# Genome ID for labeling
genome_id <- "GCF_000005845.2"

# Remove genes with ambiguous nucleotides
genes <- genes[!grepl("[^ATGCatgc]", as.character(genes))]

# Compute the length of each gene (in nucleotides)
gene_lengths <- width(genes)

# Create a codon usage table for all genes
codon_counts <- codonTable(genes)

# Calculate ENC for each gene
enc_values <- ENC(codon_counts)
head(enc_values)

# Calculate ENC prime for each gene (corrected function name)
enc_prime_values <- ENCprime(codon_counts)
# Verify that enc_prime_values is now a numeric vector
head(enc_prime_values)


# Extract the first column from the matrix
enc_prime_values <- enc_prime_values[, 1]

# Verify the extracted ENC prime values
head(enc_prime_values)



# Create a data frame with the results
genes_df <- data.frame(
  GeneID = names(genes),
  GeneLength = gene_lengths,
  ENC = enc_values,
  ENC_Prime = enc_prime_values
)

# Verify that ENC_Prime column exists
colnames(genes_df)

# Remove rows with NA values in ENC or ENC_Prime
genes_df <- genes_df %>%
  filter(!is.na(ENC), !is.na(ENC_Prime))

# Plot ENC vs Gene Length
enc_plot <- ggplot(genes_df, aes(x = GeneLength, y = ENC)) +
  geom_point(alpha = 0.3, size = 0.5, color = "blue") +
  labs(title = paste("ENC vs Gene Length for Individual Genes"),
       x = "Gene Length (bp)",
       y = "ENC") +
  theme_minimal()

# Display the plot
print(enc_plot)


# Plot ENC Prime vs Gene Length
enc_prime_plot <- ggplot(genes_df, aes(x = GeneLength, y = ENC_Prime)) +
  geom_point(alpha = 0.3, size = 0.5, color = "green") +
  labs(title = paste("ENC Prime vs Gene Length for Individual Genes"),
       x = "Gene Length (bp)",
       y = "ENC Prime") +
  theme_minimal()

# Display the plot
print(enc_prime_plot)


