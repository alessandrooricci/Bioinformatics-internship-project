#=====================================================================
# Generating boxplots for each step and put them in a grid with cowplot
#=====================================================================

# Barplot of the number of genes before and after NA filtering

before_filtering <- read_tsv("./data/dataset_cml_UnID.tsv")
before_filtering <- before_filtering %>%
  relocate(UNIPROT)

after_filtering <- read_tsv("./data/dataset_cml_filtered.tsv")

names <- c("Before_filtering", "After_filtering")
values <- c(length(before_filtering$Genes), length(after_filtering$Genes))

barplot_df <- data.frame(names, values)
barplot(height = barplot_df$values, names = barplot_df$names)

# Non-logarithmized data boxplot

not_log <- read_tsv("./data/not_log_dataframe.tsv")
boxplot(not_log)

# Logarithmized data boxplot

log_data <- read_tsv("./data/log_dataframe.tsv")
boxplot(log_data)

# Normalized data boxplot

norm_data <- read_tsv("./data/dataset_normalized.tsv")

norm_data <- norm_data %>%
  dplyr::select(-UNIPROT,-Genes)

boxplot(norm_data)

# Imputed data boxplot

imp_data <- read_tsv("./data/dataset_cml_imputed.tsv")

imp_data <- imp_data %>%
  dplyr::select(-UNIPROT,-Genes)

boxplot(imp_data)

# Put boxplots into a grid using the 'cowplot' library












