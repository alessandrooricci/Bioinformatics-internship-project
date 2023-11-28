source("0.libraries.R")
setwd("C:/Users/ricci/OneDrive/Desktop/Alessandro_tirocinio")
dataset_cml <- read_tsv("data/report.unique_genes_matrix.tsv")

#Change columns names
name_colums <- c()

for (x in c("R","S")){
   for (y in c("Ctrl","Ima")){
     for(z in c(1,2,3)) {
       nome <- paste0(x,"_",y,"_",z)
       name_colums <- c(name_colums, nome)
     }
   }
}
colnames(dataset_cml) <- name_colums
print(name_colums)
#


#
x_values <- c("R", "S")
y_values <- c("Ctrl", "Ima")
z_values <- c(1, 2, 3)

combinations <- expand.grid(X = x_values, Y = y_values, Z = z_values, stringsAsFactors = F)

for (i in 1:nrow(combinations)){
  print(paste0(combinations[i,"X"], "_", combinations[i,"Y"], "_", combinations[i,"Z"]))
}
combinations

for (i in 1:nrow(combinations)){
  if (i %% 2 == 0)
    print(paste0(combinations[i,"X"], "_", combinations[i,"Y"], "_", combinations[i,"Z"]))
}
name_columns <- apply(combinations, 1, function(row) paste0(row, collapse = "_"))
print(name_columns)
#


#Alternative way to change column names
name_colums
dataset_cml <- column_to_rownames(dataset_cml, var = "Genes")
colnames(dataset_cml) <- name_colums
dataset_cml <- rownames_to_column(dataset_cml, var = "Genes")
#













