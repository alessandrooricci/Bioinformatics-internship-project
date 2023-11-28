
# The aim of this script is to map gene symbol to UNIPROT IDs
# with three methods: UNIPROT database (ID MAPPING), clusterProfiler and 
# AnnotationDBI


#Read proteomics output of A.rename_columns.R script
dataset_cml <- read_tsv("data/proteomics_clean.tsv")

# ===================================================
# Method1: ID_mapping
# ===================================================

# Create the query to submit to UNIPROT database
gene_names <- paste(dataset_cml$Genes,collapse = " ")
gene_names

# Read the output of query on UNIPROT db
ID_mapping <- read_tsv("data/idmapping_2023_11_21.tsv")
ID_mapping <- ID_mapping %>%
  dplyr::select(-`Entry Name`,-`Protein names`)
colnames(ID_mapping) <- c("SYMBOL","UNIPROT","GENES_NAMES")

# Assign to each gene symbol multiple UNIPROT IDs on the same row
ID_mapping <- ID_mapping %>%
  group_by(SYMBOL) %>%
  reframe(UNIPROT = paste0(UNIPROT, collapse = " ; "))
ID_mapping

# 
check_3 <- inner_join(dataset_cml, ID_mapping, by = c("Genes" = "SYMBOL"))
lost_3 <- anti_join(dataset_cml, check_3, by = "Genes")
lost_3_c <- c(lost_3$Genes)
lost_3_c

#From gene names to UNIPROT ID using 
#clusterProfiler and AnnotationDbi

genes_names <- c(dataset_cml$Genes)
genes_names
BiocManager::install("org.Hs.eg.db")
BiocManager::install("AnnotationDbi")

# ===================================================
#Method2: clusterProfiler
# ===================================================

#clusterProfiler
Protein_IDs <- bitr(genes_names, fromType = "SYMBOL", toType = "UNIPROT", OrgDb = "org.Hs.eg.db")

cluster_prof <- Protein_IDs %>%
  group_by(SYMBOL) %>%
  reframe(UNIPROT = paste0(UNIPROT, collapse = " ; "))
cluster_prof

check_1 <- inner_join(dataset_cml, cluster_prof, by = c("Genes" = "SYMBOL"))
lost_1 <- anti_join(dataset_cml, check_1, by = "Genes")
lost_1_c <- c(lost_1$Genes)
lost_1_c

# ===================================================
#Method3: AnnotationDbi
# ===================================================


#AnnotationDbi
uniKeys <- keys(org.Hs.eg.db, keytype="UNIPROT")
gene_ann <- select(org.Hs.eg.db, keys = genes_names, keytype = "SYMBOL", columns = c("SYMBOL", "UNIPROT"))

ann_Dbi <- Protein_IDs %>%
  group_by(SYMBOL) %>%
  reframe(UNIPROT = paste0(UNIPROT, collapse = " ; "))
ann_Dbi

check_2 <- inner_join(dataset_cml, ann_Dbi, by = c("Genes" = "SYMBOL"))
lost_2 <- anti_join(dataset_cml, check_2, by = "Genes")
lost_2_c <- c(lost_2$Genes)
lost_2_c

setdiff(lost_1_c, lost_2_c) # Annotation DBI and ClusterProfiler output are the same
# The best is UNIPROT ID mapping

# ===================================================
# Correct UNIPROT ID mapping
# ===================================================

# Transform manually wrong gene symbols
#Deleting lost_3 genes from dataset_cml
dataset_cml_1 <- dataset_cml[!(dataset_cml$Genes %in% lost_3_c), ]

#Changing rows name in lost_3
lost_3$Genes <- c(paste0("MARCHF", c(2,5,6,7)), "NCBP2AS2","PPP5D1P","SEP15", paste0("SEPT", c(1,10,11,2,5,6,7,8,9)))

#Join lost_3 to dataset_cml 
dataset_cml_clean <- bind_rows(dataset_cml_1, lost_3)

# Create the new query
lost_3_gene <- paste(lost_3$Genes,collapse = " ")
lost_3_gene

# Read the output of query on UNIPROT db
ID_mapping_lost <- read_tsv("data/idmapping_2023_11_22.tsv")
ID_mapping_lost <- ID_mapping_lost %>%
  dplyr::select(-`Entry Name`,-`Protein names`)
colnames(ID_mapping_lost) <- c("SYMBOL","UNIPROT","GENES_NAMES")

#Delete PPP5D1 gene from main dataset
dataset_cml_clean <- dataset_cml_clean[dataset_cml_clean$Genes != "PPP5D1P", ]

#Unify ID_mapping datasets with all UNIPROT IDs
ID_mapping <- bind_rows(ID_mapping, ID_mapping_lost)

#Add UNIPROT ID column to clean dataset cml
dataset_cml_ID <- cbind(dataset_cml_clean, UNIPROT = ID_mapping$UNIPROT)

# Save a final df
write_tsv(dataset_cml_ID, "./data/dataset_cml_UnID.tsv")

