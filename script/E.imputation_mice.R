# MICE (Multivariate Imputation by Chained Equations) is used 
# in this script to deal with missing data

source("./script/0.libraries.R")
proteomics_norm <- read_tsv("./data/dataset_normalized.tsv")

# Remove columns that are not needed
proteomics_imp <- proteomics_norm %>%
  dplyr::select(-UNIPROT, -Genes)

# Init MICE
init = mice(proteomics_imp, maxit = 0)
meth = init$method
predM = init$predictorMatrix

# Choose method 
meth[1:length(meth)] <- "norm"

# Run the multiple (m=5) imputation
imputed = mice(proteomics_imp, method = meth , predictorMatrix = predM, m = 5)

# Create a dataset after imputation
data_imputed <- complete(imputed)

# Add deleted columns
data_imputed <- bind_cols(data_imputed, proteomics_norm %>% 
                            dplyr::select(UNIPROT, Genes)) %>%
  relocate(UNIPROT, Genes) %>%
  as_tibble()

# Save imputed data frame
write_tsv(data_imputed, "./data/dataset_cml_imputed.tsv")

# Imputed data boxplot
boxplot_imp <- data_imputed %>%
  dplyr::select(-UNIPROT, -Genes)
boxplot(boxplot_imp) 
