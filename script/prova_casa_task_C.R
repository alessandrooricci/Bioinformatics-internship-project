
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
## Create tibbles and change names

NA_values_per_group_1 <- as.tibble(rowSums(is.na(proteomics_excercises[, 2:4])))
colnames(NA_values_per_group_1) <- c("NA_group_1")

NA_values_per_group_2 <- as.tibble(rowSums(is.na(proteomics_excercises[, 5:7])))
colnames(NA_values_per_group_2) <- c("NA_group_2")

NA_values_per_group_3 <- as.tibble(rowSums(is.na(proteomics_excercises[, 8:10])))
colnames(NA_values_per_group_3) <- c("NA_group_3")

NA_values_per_group_4 <- as.tibble(rowSums(is.na(proteomics_excercises[, 11:13])))
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
  filter(rowSums(NA_matrix) <= 10)
# Opposite to line 44
NA_table_filtered <- rownames_to_column(NA_table_filtered, "Genes")

#====================================================================================
# PROVA
# Create 4 groups to filter using flag columns applying the condition

R_Ctrl <- proteomics_excercises[,3:5]
R_Ctrl$flag1 <- rowSums(!is.na(R_Ctrl)) >= 2

R_Ima <- proteomics_excercises[,6:8]
R_Ima$flag2 <- rowSums(!is.na(R_Ima)) >= 2

S_Ctrl <- proteomics_excercises[,9:11]
S_Ctrl$flag3 <- rowSums(!is.na(S_Ctrl)) >= 2

S_Ima <- proteomics_excercises[,12:14]
S_Ima$flag4 <- rowSums(!is.na(S_Ima)) >= 2

output <- cbind(proteomics_excercises[,1:2],R_Ctrl,R_Ima,S_Ctrl,S_Ima)

output <- output[rowSums(output[, c("flag1", "flag2", "flag3", "flag4")]) >= 1, ]

output <- output %>%
  dplyr::select(-flag1, -flag2, -flag3, -flag4)

#====================================================================================

# Turn 'Genes' row to column 
NA_table_filtered <- rownames_to_column(NA_table_filtered, "Genes")

#=======================================
# Join gene names to the main data set 
#=======================================

proteomics_filtered <- semi_join(proteomics_excercises,NA_table_filtered , by = "Genes")

# Save the final data set filtered
write_tsv(proteomics_filtered, "./data/dataset_cml_filtered.tsv")