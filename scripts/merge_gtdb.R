#MERGE 180 gtdb220

#Load Required Packages
library(dplyr)

#Load the Data Frames
# Load codon statistics (df_stats) and growth data (growth_df)
df_stats <- readRDS("genome_predictions_truncated.rds")
load("/gpfs/projects/WeissmanGroup/gtdb220/gtdb220_genome_predictions_full.rda")  # file contains an object named 'growth_df'


#only include the columns Genome and doubling time (d):
pred_df <- pred_df %>% select(Genome, d)
pred_df <- pred_df %>% mutate(Genome = gsub("\r", "", Genome))


#clean the Genome Column in df_stats
df_stats <- df_stats %>%
  select(Genome, CUBHE, ConsistencyHE, CPB)

#remove the suffix from the Genome names
df_stats <- df_stats %>%
  mutate(
    Genome = as.character(Genome),
    #trim whitespace (spaces, tabs, newlines) from both ends
    Genome = trimws(Genome),
    #remove carriage returns if present
    Genome = gsub("\r", "", Genome),
    #remove the exact suffix "_genomic.ffn"
    Genome = gsub("_genomic\\.ffn$", "", Genome)
  )


#merge with growth data
#keep only the Genome and d columns from pred_df
merged_df <- merge(
  df_stats,
  pred_df %>% select(Genome, d),
  by = "Genome",
  all.x = TRUE
)
merged_df <- merged_df[complete.cases(merged_df), ]


cat("Merged data frame has", nrow(merged_df), "rows.\n")
cat("Columns in merged data:\n")
print(names(merged_df))
print(head(merged_df))

#save the Merged Data Frame
saveRDS(merged_df, file = "merged_codon_stats_growth.rds")
write.csv(merged_df, file = "merged_codon_stats_growth.csv", row.names = FALSE)

cat("Merged data frame saved as 'merged_codon_stats_growth.rds' and 'merged_codon_stats_growth.csv'.\n")

