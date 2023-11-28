source("./script/0.libraries.R")

#==========================
# Excercises
#==========================

# Load proteomics df
proteomics_excercises <- read_tsv("./data/proteomics_clean.tsv")

# Find total NA values in df
NA_values <- sum(is.na(proteomics_excercises))

# NA in each column
NA_values_col <- colSums(is.na(proteomics_excercises))

# NA in each row
NA_values_row <- rowSums(is.na(proteomics_excercises))

# NA values total counting
sum(is.na(proteomics_excercises))

# Rows with no NA values
clean_rows <- proteomics_excercises[complete.cases(proteomics_excercises), ]

# Rows with one or none NA values
one_or_none_rows <- proteomics_excercises %>%
  filter(rowSums(is.na(proteomics_excercises)) <= 1)

#===============================================================================

# Find NA values in each group/condition
## Creating a tibble with gene names

gene_names <- as.tibble(proteomics_excercises$Genes)

colnames(gene_names) <- c("Genes") # Change name

## Join tibbles with their NA values
### Create tibbles and change names

NA_values_per_group_1 <- as.tibble(rowSums(is.na(proteomics_excercises[, 2:4])))
colnames(NA_values_per_group_1) <- c("NA_group_1")

NA_values_per_group_2 <- as.tibble(rowSums(is.na(proteomics_excercises[, 5:7])))
colnames(NA_values_per_group_2) <- c("NA_group_2")

NA_values_per_group_3 <- as.tibble(rowSums(is.na(proteomics_excercises[, 8:10])))
colnames(NA_values_per_group_3) <- c("NA_group_3")

NA_values_per_group_4 <- as.tibble(rowSums(is.na(proteomics_excercises[, 11:13])))
colnames(NA_values_per_group_4) <- c("NA_group_4")

### Join tibbles to create a unique table

NA_table <- cbind(gene_names, NA_values_per_group_1, NA_values_per_group_2, 
                  NA_values_per_group_3, NA_values_per_group_4)

#===================================================================================
# Filter the table by keeping genes that have at least in a group 2 out of 3 values
#===================================================================================

NA_matrix <- NA_table %>%
  column_to_rownames("Genes")

NA_table_filtered <- NA_table %>%
  filter(rowSums(NA_matrix) <= 10)

# devo prendere i geni che hanno non più di 10 NA values 



#PROVA

prova <- as.tibble(rowSums(is.na(str_detect(proteomics_excercises, pattern = "R_Ctrl"))))
colnames(prova) <- c("NA_group_1")


