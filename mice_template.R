
# MICE code

# Read the dataframe to impute
phospho_df <- read_tsv('./processed_input/cleaning/phosphoproteome_clean.tsv')

# Remove column that aren't samples
phospho_df_imp <- phospho_df %>% dplyr::select(-c(UNIPROT, gene_name, position, aminoacid, phosphopeptide, num_peps))

# If column names contains '-' it could give an error, so be careful
colnames(phospho_df_imp) <- c(paste0('Patient', 1:47))

# Init mice
init = mice(phospho_df_imp, maxit=0)
meth = init$method
predM = init$predictorMatrix

# Choose the method to use
meth[1:length(meth)] <- "norm"

# Run the multiple (m=5) imputation.
imputed = mice(phospho_df_imp, method= meth , predictorMatrix=predM, m=5)

#Create a dataset after imputation.
imputed <- complete(imputed)

# Add the removed columns
phospho_imputed_df <- bind_cols( imputed,
                                 phospho_df %>% dplyr::select(UNIPROT, gene_name, position, aminoacid, phosphopeptide, num_peps)) %>%
  relocate(UNIPROT, gene_name, position, aminoacid, phosphopeptide, num_peps) %>%
  as_tibble()

# If you change the column names give the old names again
colnames(phospho_imputed_df) <- colnames(phospho_df)

# Save the output.
write_tsv(phospho_imputed_df, './processed_input/imputation/phospho_imputed_df.tsv')
