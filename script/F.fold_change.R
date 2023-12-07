# In this script we want to calculate the fold change between 
# R_Ctrl cells and S_Ctrl cells and create 
# a density graph by identifying the median

# Load libraries
source("./script/0.libraries.R")

# Load data frame
data_imputed <- read_tsv("./data/dataset_cml_imputed.tsv")

# Divide data frame into the two groups of interest
R_Ctrl_1 <- data_imputed %>%
  dplyr::select(starts_with("R_Ctrl"))

S_Ctrl_1 <- data_imputed %>%
  dplyr::select(starts_with("S_Ctrl"))

#Calculate fold change by subtracting the average of the rows of the two data frames
fold_change <- as.tibble(rowMeans(R_Ctrl_1) - rowMeans(S_Ctrl_1))

# Change column name
fold_change <- fold_change %>%
  dplyr::rename(Fold_Change_values = value)  

# Bind fold change column to main data frame
data_imputed <- cbind(data_imputed, fold_change)

#=============================================
# Density plot 
#=============================================

#Calculate median value
median_value <- median(fold_change$Fold_Change_values)
# Density plot building
fold_change %>%
  ggplot( aes(x = Fold_Change_values)) +
  geom_density(fill = "#30D5C8", alpha = 0.8) +
  geom_vline(xintercept = median_value, color = "#FF0000")

#=============================================

# Find the number of genes that have a fold-change in absolute value > 2
genes_greater <- data_imputed %>%
  dplyr::filter(abs(data_imputed$Fold_Change_values) > 2)






