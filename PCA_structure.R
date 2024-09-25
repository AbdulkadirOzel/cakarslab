library(limma)       # Load limma for linear modeling and differential expression analysis
library(Biobase)     # Load Biobase for handling biological data
library(ggplot2)     # Load ggplot2 for data visualization
library(cluster)     # Load cluster for clustering algorithms
library(factoextra)  # Load factoextra for visualization of clustering results
library(ggfortify)   # Load ggfortify for PCA visualization
library(cowplot)     # Load cowplot for combining plots


setwd("give_path")   # Set working directory (provide the appropriate path)

# Load data
SDRF <- read.delim("give_path/samples.txt", check.names=FALSE, stringsAsFactors=FALSE) 
# Read sample metadata from a file

x <- read.maimages(SDRF[,"Source_Name"], source="agilent", green.only=TRUE, other.columns="gIsWellAboveBG") 
# Read microarray data from Agilent arrays (green channel only)

y <- normalizeBetweenArrays(x, method="quantile") 
# Normalize the microarray data between arrays using quantile normalization

# Create data frame
a <- as.data.frame(y$E) 
# Convert the expression matrix (normalized intensities) to a data frame

new_column_names <- sub(".*/", "", colnames(a)) 
# Simplify the column names by removing directory paths, leaving only file names

colnames(a) <- new_column_names 
# Apply the simplified column names to the data frame

# PCA
pca_result <- prcomp(t(a)) 
# Perform Principal Component Analysis (PCA) on the transposed data (samples as rows)

# Check the column names
cat("Column names:\n") 
# Display the column names for reference

print(colnames(a)) 
# Print the column names to verify

# Identify reference strain data
reference_pattern <- "WT"  # Adjust the pattern to match the naming of reference data
reference_indices <- grep(reference_pattern, colnames(a))  
# Identify the indices of reference strains in the data based on the pattern

reference_strains <- pca_result$x[reference_indices, ]  
# Extract PCA results for the reference strains

remaining_data <- pca_result$x[-reference_indices, ]  
# Extract PCA results for the non-reference (remaining) strains

# Check dimensions
cat("Number of reference strains: ", length(reference_indices), "\n") 
# Print the number of reference strains

cat("Dimensions of remaining data: ", dim(remaining_data), "\n") 
# Print the dimensions of the remaining data

cat("Dimensions of reference strains: ", dim(reference_strains), "\n") 
# Print the dimensions of the reference strains

# K-means clustering - Exclude reference strain data
set.seed(123) 
# Set random seed for reproducibility

num_clusters <- 5 
# Define the number of clusters for k-means clustering

if (nrow(remaining_data) >= num_clusters) {
  kmeans_result <- kmeans(remaining_data[, 1:2], centers=num_clusters)  
  # Apply k-means clustering on the first two principal components for the remaining data
  print(kmeans_result) 
  # Print the k-means result
} else {
  cat("Insufficient data points for k-means clustering with", num_clusters, "clusters.\n") 
  # Display a message if there are not enough data points for clustering
}

remaining_clusters <- kmeans_result$cluster 
# Extract cluster assignments for the remaining data

# Group information
group_colors <- c("red", "green", "blue", "orange", "pink", "black", "cyan", "purple", "gold", "brown", "gray", "magenta", "chartreuse4") 
# Define colors for each group

group_names <- c("Caffeine-resistant strain", "Cobalt-resistant strain", "Coniferylaldehyde-resistant strain", "Ethanol-resistant (B2) strain", 
                 "Ethanol-resistant (B8) strain", "Reference strain", "Iron-resistant strain", "Long-lived strain", "Nickel-resistant strain", 
                 "Oxidative stress-resistant strain", "2-phenylethanol-resistant strain", "Silver-resistant strain", "Boron-resistant strain") 
# Define descriptive group names

group_indices <- c("caf", "cobalt", "con-al", "b2", "b8", "Ref", "iron", "LL", "nic", "oxi", "pe", "sil", "bor") 
# Define patterns to match group names to sample names

# Create a data frame for ggplot
remaining_pca_data <- data.frame(PC1 = remaining_data[, 1], PC2 = remaining_data[, 2], Cluster = remaining_clusters, Sample = colnames(a)[-reference_indices]) 
# Create a data frame with the first two principal components, cluster assignments, and sample names for the remaining strains

remaining_pca_data$Group <- NA  
# Initialize the 'Group' column as NA

for (i in 1:length(group_indices)) {
  match_indices <- grep(paste0("^", group_indices[i]), remaining_pca_data$Sample)  
  # Find samples matching the group indices
  remaining_pca_data$Group[match_indices] <- group_names[i]  
  # Assign group names to the matching samples
}

# Handle NA values in the data
remaining_pca_data <- na.omit(remaining_pca_data)  
# Remove rows with NA values

# Reorder cluster labels based on mean PC1 values
cluster_centroids <- aggregate(PC1 ~ Cluster, data = remaining_pca_data, FUN = mean)  
# Calculate the mean PC1 values for each cluster

ordered_clusters <- cluster_centroids[order(cluster_centroids$PC1), "Cluster"]  
# Order the clusters based on their mean PC1 values

# Create a mapping from old cluster numbers to new ordered cluster numbers
cluster_mapping <- setNames(seq_along(ordered_clusters), ordered_clusters) 
# Map old cluster numbers to new ordered cluster numbers

# Apply new cluster numbers to the data
remaining_pca_data$OrderedCluster <- factor(remaining_pca_data$Cluster, levels = ordered_clusters, labels = seq_along(ordered_clusters)) 
# Apply the new ordered cluster numbers

# Define custom shapes for clusters
cluster_shapes <- c(0, 1, 2, 3, 4, 5, 6)  
# Define shapes for each cluster

# Combine reference strain data
reference_pca_data <- data.frame(PC1 = reference_strains[, 1], PC2 = reference_strains[, 2], Cluster = "Reference", Sample = colnames(a)[reference_indices])  
# Create a data frame for the reference strains

reference_pca_data$Group <- "Reference strain"  
# Assign 'Reference strain' as the group name for reference strains

# Ensure columns match for combining data frames
reference_pca_data$OrderedCluster <- "Reference"  
# Set the 'OrderedCluster' column for reference strains

reference_pca_data <- reference_pca_data[, colnames(remaining_pca_data)]  
# Reorder columns to match the remaining PCA data

# Combine data for plotting
combined_pca_data <- rbind(remaining_pca_data, reference_pca_data)  
# Combine the remaining PCA data with the reference strain data

# Elbow method plot
elbow_plot <- fviz_nbclust(remaining_data[, 1:2], kmeans, method = "wss") + 
  geom_vline(xintercept = 5, linetype = 2) +  
  theme_minimal() + 
  theme(
    plot.margin = margin(2, 2, -1, 1, "cm"), 
    axis.title = element_text(size = 8), 
    axis.text = element_text(size = 8), 
    title = element_text(size = 0)
  )  
# Create an elbow plot to determine the optimal number of clusters and add a vertical line at 5 clusters

# PCA plot with clusters
p <- ggplot(combined_pca_data, aes(x = PC1, y = PC2, color = Group, shape = Cluster)) +
  geom_point(size = 3, stroke = 1.5) +  # Increase point border thickness
  scale_color_manual(values = c(group_colors, "black")) +  # Include color for reference strain
  scale_shape_manual(values = c(cluster_shapes, 8)) +  # Include shape for reference strain
  labs(x = "PC1 (46.54%)", y = "PC2 (14.67%)", shape = "Cluster", color = "Group") +
  theme_bw() + 
  theme(legend.position = "right")  
# Create a PCA plot showing the clusters and reference strains with custom colors and shapes

# Combine the two plots
combined_plot <- ggdraw() +
  draw_plot(p) +
  draw_plot(elbow_plot, x = 0.02, y = 0.75, width = 0.35, height = 0.35)  
# Combine the PCA plot with the elbow plot

# Display the combined plot
print(combined_plot)  
# Show the final combined plot
