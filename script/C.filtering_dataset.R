
## In this script we want to filter our data set by removing genes 
## that do not have at least in a group 2/3 values

source("./script/0.libraries.R")
proteomics_excercises <- read_tsv("./data/dataset_cml_UnID.tsv") %>%
  relocate("UNIPROT")

#===================================================================================
# Dataset's filtering by keeping genes that have at least in a group 2/3 values
#===================================================================================

# Find NA values in each group/condition
## Creating two different tibbles: one with only gene names and another 
## one with gene names and UNIPROT ID

gene_names_UNID <- data.frame(Genes = proteomics_excercises$Genes,
                              UNIPROT = proteomics_excercises$UNIPROT)

only_gene <- as.tibble(proteomics_excercises$Genes)
colnames(only_gene) <- c("Genes")

# Join tibbles with their NA values
## Create tibbles, search by pattern and change columns names

group1_id <- grepl('R_Ctrl', colnames(proteomics_excercises))
NA_values_per_group_1 <- as.tibble(rowSums(is.na(proteomics_excercises[, group1_id])))
colnames(NA_values_per_group_1) <- c("NA_group_1")

group2_id <- grepl('R_Ima', colnames(proteomics_excercises))
NA_values_per_group_2 <- as.tibble(rowSums(is.na(proteomics_excercises[, group2_id])))
colnames(NA_values_per_group_2) <- c("NA_group_2")

group3_id <- grepl('S_Ctrl', colnames(proteomics_excercises))
NA_values_per_group_3 <- as.tibble(rowSums(is.na(proteomics_excercises[, group3_id])))
colnames(NA_values_per_group_3) <- c("NA_group_3")

group4_id <- grepl('S_Ima', colnames(proteomics_excercises))
NA_values_per_group_4 <- as.tibble(rowSums(is.na(proteomics_excercises[, group4_id])))
colnames(NA_values_per_group_4) <- c("NA_group_4")

# Join tibbles to create a unique table
NA_table <- cbind(only_gene, NA_values_per_group_1, NA_values_per_group_2, 
                  NA_values_per_group_3, NA_values_per_group_4)

#=============================
# Dataset's filtering
#=============================

# Turn 'Genes' column to row and create a new variable
NA_matrix <- NA_table %>%
  column_to_rownames("Genes")
# Using filter() to filter genes of interest  
NA_table_filtered <- NA_matrix %>%
  filter(rowSums(NA_matrix <= 1) >= 1) 


# Turn 'Genes' row to column
NA_table_filtered <- rownames_to_column(NA_table_filtered, "Genes")


#=======================================
# Join gene names to the main data set 
#=======================================

proteomics_filtered <- semi_join(proteomics_excercises,NA_table_filtered , by = "Genes")

# Save the final data set filtered
write_tsv(proteomics_filtered, "./data/dataset_cml_filtered.tsv")

