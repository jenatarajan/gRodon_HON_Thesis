#This script is for cutting genes to 30 base pairs, the process was repeated for 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 390, 420, 450, 480, 510 bp
#loading libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("Biostrings", quietly = TRUE))
  BiocManager::install("Biostrings")
if (!requireNamespace("parallel", quietly = TRUE))
  install("parallel")
library(Biostrings)
library(parallel)
#paths of input and output folder
input_folder <- "/gpfs/projects/WeissmanGroup/gRodon_training/madin_genomes"
output_folder <- "/gpfs/projects/WeissmanGroup/gRodon_training/all_short_madin_genomes_30"
#create output folder
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}
#list of all genome files in the input folder
genome_files <- list.files(input_folder, full.names = TRUE, pattern = "\\.(fasta|fa|fna|gz)$")
#function to cut the sequences to 30 base pairs
geneSubseq <- function(my_seq){
  if (nchar(my_seq) >= 30) {
    subseq(my_seq, start = 1, end = 30)
  } else {
    NA  #skip sequences shorter
  }
}
genomeTrim <- function(genome_file,output_folder){
  genome_sequences <- readDNAStringSet(genome_file)
  #truncate each sequence to the first base pairs
  truncated_sequences <- lapply(genome_sequences, geneSubseq)
  truncated_sequences <- DNAStringSet(truncated_sequences[!is.na(truncated_sequences)])
  #output file path
  output_file <- file.path(output_folder, paste0("truncated_", basename(genome_file)))
  #write truncated sequences to output file path
  writeXStringSet(truncated_sequences, filepath = output_file, format = "fasta")
  print(paste("Truncated", basename(genome_file), "to 30 bp and saved to", output_file))
}
#Run in parallel on many genomes
mclapply(X = genome_files,
         FUN = genomeTrim,
         output_folder = output_folder,
         mc.cores = 28) 
