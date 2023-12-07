# DEP code
# This package provides an integrated analysis workflow for 
# robust and reproducible analysis of mass spectrometry 
# proteomics data for differential protein expression or differential enrichment. 

# Load data 
init_data <- read_tsv("./data/dataset_cml_UnID.tsv")
# Select columns 
init_data <- init_data %>%
  relocate(UNIPROT) %>%
  dplyr::select(-starts_with('R_Ima'), -starts_with('S_Ima'))

# Generate data with 'name' and 'ID' columns 
data_unique <- make_unique(init_data, "Genes", "UNIPROT", delim = ";")
# Personalization
data_unique$UNIPROT <- NULL
data_unique$Genes <- NULL
data_unique <- data_unique %>%
  relocate(name, ID)

# Generate a SummarizedExperiment object using an experimental design
proteomics_expdesign <- data.frame(label = c(paste0("R_Ctrl_", 1:3), 
                                             paste0("S_Ctrl_", 1:3)),
                                   condition = c(rep("R_Ctrl", each = 3), 
                                                 rep("S_Ctrl", each = 3)),
                                   replicate = c(1:3,1:3))

R_Ctrl_columns <- grep("R_Ctrl", colnames(init_data)) # get R_Ctrl columns numbers
S_Ctrl_columns <- grep("S_Ctrl", colnames(init_data)) # get S_Ctrl columns numbers
columns <- c(R_Ctrl_columns, S_Ctrl_columns)
experimental_design <- proteomics_expdesign
data_se <- make_se(data_unique, columns, experimental_design) # Create the SE
data_se 

# Plot a barplot of the protein identification overlap between samples
plot_frequency(data_se) 

# Filter for proteins that are identified 
# in 2 out of 3 replicates of at least one condition
data_filt <- filter_missval(data_se, thr = 1)
plot_frequency(data_filt)

# Plot a barplot of the number of identified proteins per samples
plot_numbers(data_filt)

# Plot a barplot of the protein identification overlap between samples
plot_coverage(data_filt)

#===============================================================================
# DATA NORMALIZING
#===============================================================================

data_norm <- normalize_vsn(data_filt)

# Visualize normalization by boxplots for all samples 
# before and after normalization
plot_normalization(data_norm, data_filt)

#===============================================================================
# DATA IMPUTING
#=============================================================================== 

# Plot a heatmap of proteins with missing values
plot_missval(data_filt)

# Plot intensity distributions and cumulative fraction 
# of proteins with and without missing values
plot_detect(data_filt)

#cowplot::plot_grid(plotlist = list(a1,a2), ncol = 2)

# Impute missing data using random draws from a Gaussian 
# distribution centered around a minimal value 
# (for MNAR; missing not a random)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
plot_imputation(data_norm, data_imp)

# Differential enrichment analysis  based on linear models and empherical Bayes statistics

# Test every sample versus control
data_diff <- test_diff(data_imp, type = "control", control = "S_Ctrl")

# Denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(2))


# Plot the first and second principal components (?)
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)

# Plot the Pearson correlation matrix
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")

# Plot a heatmap of all significant proteins with the data centered per protein
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = FALSE,
             indicate = c("condition", "replicate"))

# Plot a heatmap of all significant proteins (rows) and the tested contrasts (columns)
plot_heatmap(dep, type = "contrast", kmeans = TRUE, 
             k = 6, col_limit = 10, show_row_names = FALSE)

# Plot a volcano plot for the contrast "R_Ctrl vs S_Ctrl""
plot_volcano(dep, contrast = "R_Ctrl_vs_S_Ctrl", label_size = 2, add_names = T)









