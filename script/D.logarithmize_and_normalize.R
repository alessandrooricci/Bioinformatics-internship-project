# In this script we want to change the values in the table 
# and put them in logarithmic form (log2 x + 1), then we impute with 'mice' 

source("./script/0.libraries.R")
proteomics <- read_tsv("./data/dataset_cml_filtered.tsv")

#===========================================================
# Data logarithmizing and data normalizing
# ==========================================================

normalized_data <- proteomics %>%
  # Select only the columns of interest (excluding 'UNIPROT' and 'GENES')
  select_if(is.numeric) %>%
  # Apply logarithmic transformation
  mutate(across(everything(), ~log2(. + 1))) %>%
  # Apply normalization dividing by the median
  mutate(across(everything(), ~./median(., na.rm = TRUE)))

# Join the 'UNIPROT' and 'GENES' columns to the normalized dataframe
data_normalized <- cbind(proteomics[, c("UNIPROT", "Genes")], normalized_data)

# Creating data frame with logaritmic values for boxplot
logaritmic_data <- proteomics %>%
  dplyr::select(-UNIPROT, -Genes) %>%
  mutate(across(everything(), ~log2(. + 1)))

# Box plot with logaritmic values
box_log <- boxplot(logaritmic_data)
# BOx plot with normalized values
box_norm <- boxplot(normalized_data)

write_tsv(data_normalized, "./data/dataset_normalized.tsv")



