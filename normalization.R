library(limma)    # Load the limma package for linear modeling and differential expression analysis
library(Biobase)  # Load Biobase package, which is necessary for working with expression data
library(GEOquery) # Load GEOquery for retrieving data from the Gene Expression Omnibus (GEO)

setwd("give_path") # Set working directory (path should be provided)

# Limma user guide
# Inspired by section 17.4 "Time Course Effects of Corn Oil on Rat Thymus with Agilent 4x44K Arrays"

SDRF <- read.delim("give_path/E-GEOD-33005.sdrf.txt", check.names=FALSE, stringsAsFactors=FALSE) 
# Load sample description file (SDRF) to obtain metadata
# It contains file names and their status like wild-type or a strain

design_1 <- model.matrix(~0 + Status, data=SDRF) 
# Create a design matrix for linear modeling, based on the 'Status' column in SDRF

colnames(design_1) <- c('mutant', 'ref') 
# Assign column names to the design matrix for 'mutant' and 'ref' groups

rownames(design_1) <- c("caf1", "caf2", "caf3", "W1", "W2", "W3") 
# Assign row names (sample identifiers) to the design matrix
# This is for caffeine-resistant strain. Of course, the designs can be changed.

cm <- makeContrasts(status = mutant - ref, levels = design_1) 
# Create contrast matrix to compare 'mutant' vs 'ref' groups

x <- read.maimages(SDRF[,"Source_Name"], 
                   source="agilent", green.only=TRUE, other.columns="gIsWellAboveBG") 
# Read microarray data from Agilent arrays; use green channel only and include signal quality information

y <- x 
# Assign raw data to variable 'y' (background correction skipped)

y <- normalizeBetweenArrays(y, method="quantile") 
# Normalize data between arrays using quantile normalization

plot(density(y$E[,1])) 
# Plot the density of the first sample after normalization

plot(density(x$E[,1])) 
# Plot the density of the first sample before normalization (raw data)

boxplot(y$E, las=2) 
# Generate boxplot to visualize distribution of expression values for each sample

dim(y) 
# Check the dimensions of the 'y' object (number of features and samples)

# Should also check with a histogram

# Control <- y$genes$ControlType==1L 
# Filter control probes (commented out)

# IsExpr <- rowSums(y$other$gIsWellAboveBG > 0) >= 4 
# Define 'expressed' probes based on signal quality (commented out)

# yfilt <- y[!Control & IsExpr, ] 
# Filter out control probes and low-expressed probes (commented out)

# dim(yfilt) 
# Check dimensions of the filtered dataset (commented out)

fit <- lmFit(y, design_1) 
# Fit a linear model to the expression data based on the design matrix

fit2 <- contrasts.fit(fit, contrasts = cm) 
# Apply the contrast matrix to the fitted model

fit2 <- eBayes(fit2, trend=TRUE, robust=TRUE) 
# Perform empirical Bayes moderation with trend and robust options for more reliable estimates

results <- decideTests(fit2) 
# Determine significantly differentially expressed genes

summary(results) 
# Summarize the results (number of upregulated/downregulated/unchanged genes)

toptable <- topTable(fit2, coef = 1, number = Inf) 
# Generate a table of the top differentially expressed genes based on the contrast

plot(toptable$logFC) 
# Plot the log-fold changes (logFC) of all genes

plot(abs(toptable$logFC)) 
# Plot the absolute values of log-fold changes

table(abs(toptable$logFC) >= 2) 
# Create a table of genes with absolute log-fold change greater than or equal to 2

aver_topTable <- avereps(toptable, ID = toptable$ProbeName) 
# Average the expression values for genes represented by multiple probes

write.csv(aver_topTable, "give_path/filename.tsv", row.names=FALSE, quote=FALSE) 
# Write the results to a CSV file (path should be provided)
