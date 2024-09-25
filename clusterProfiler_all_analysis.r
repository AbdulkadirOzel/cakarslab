
library(clusterProfiler)
library(DOSE)
library(ggplot2)
setwd("path")

# Reading in all the TSV files for various conditions
caffeine = read.csv("path/caffeine.tsv", sep=',', header=TRUE)
b8_ethanol = read.csv("path/B8-ethanol.tsv", sep=',', header=TRUE)
coniferylaldehyde = read.csv("path/coniferyl-aldehyde.tsv", sep=',', header=TRUE)
silver = read.csv("path/silver.tsv", sep=',', header=TRUE)
iron = read.csv("path/iron.tsv", sep=',', header=TRUE)
oxidative = read.csv("path/oxidative.tsv", sep=',', header=TRUE)
phenylethanol = read.csv("path/phenylethanol.tsv", sep = ',', header = TRUE)
b2_ethanol = read.csv("path/B2-ethanol.tsv", sep=',', header=TRUE) 
boron = read.csv("path/boron.tsv", sep=',', header=TRUE)
long_lived  = read.csv("path/lifespan.tsv", sep=',', header=TRUE)
cobalt = read.csv("path/cobalt.tsv", sep=',', header=TRUE)
nickel = read.csv("path/nickel.tsv", sep=',', header=TRUE)

# Store all the dataframes in a list for easy iteration
dataframes <- list(
  caffeine = caffeine,
  b8_ethanol = b8_ethanol,
  coniferylaldehyde = coniferylaldehyde,
  silver = silver,
  iron = iron,
  oxidative = oxidative,
  phenylethanol = phenylethanol,
  b2_ethanol = b2_ethanol,
  boron = boron,
  long_lived = long_lived,
  cobalt = cobalt,
  nickel = nickel
)

# Vector to collect all unique pathway names across all samples
all_pathways <- c()

# List to store NES (Normalized Enrichment Scores) values for each sample
nes_list <- list()

# Loop through each dataframe in the list
for (df_name in names(dataframes)) {
  df <- dataframes[[df_name]]
  
  # Filtering rows where P.Value is less than 0.05
  df <- subset(df, P.Value < 0.05)
  
  # Create a named gene list using logFC and set gene names as names
  genelist <- df[, "logFC"]
  names(genelist) <- as.character(df[, "SystematicName"])
  
  # Sort gene list in decreasing order of logFC
  genelist <- sort(genelist, decreasing = TRUE)
    
  # Perform KEGG pathway enrichment analysis
  kk <- gseKEGG(geneList = genelist, organism = 'sce')
  
  # Collect unique pathway names
  all_pathways <- unique(c(all_pathways, kk@result$Description))
  
  # Collect NES (Normalized Enrichment Score) for each pathway
  nes_list[[df_name]] <- kk@result %>% 
    select(Description, NES) %>% 
    rename(Pathway = Description, !!df_name := NES)
}

# Create a new dataframe to store NES values for each sample with pathways as rows
nes_df <- data.frame(Pathway = all_pathways)

# Initialize NES values for each sample with 0
for (sample in names(dataframes)) {
  nes_df[[sample]] <- 0
}

# Fill in NES values in the dataframe for each sample
for (sample in names(dataframes)) {
  nes_data <- nes_list[[sample]]
  
  # Loop through the NES data and update corresponding rows in the NES dataframe
  for (i in 1:nrow(nes_data)) {
    pathway <- nes_data$Pathway[i]
    nes_value <- nes_data[[sample]][i]
    nes_df[nes_df$Pathway == pathway, sample] <- nes_value
  }
}

# Write the final NES values matrix to a CSV file
write.csv(nes_df, "path/NES_values_matrix.csv", row.names = FALSE)