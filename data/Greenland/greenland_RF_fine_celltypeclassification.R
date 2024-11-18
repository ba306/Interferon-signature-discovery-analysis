# Cell type labeling using RF trained on our immune cell type genes
# Gene list: our immune cell type genes 
# Reference: Hao
# Query: Greenland

library(foreach)
library(qs)
library(dplyr)
set.seed(42)


# Call the function for the training and classifier
source("./functions/RF_celltype_classifier.R")
source("./functions/benchmarking_prediction_metrics.R")

#### Define labeling of the outpus ####

# How to name the results 

ref_name="Hao_fine"
query_name="Greenland"
save_dir= "./data/Greenland/"

# Colname for the celltype metadata in the reference

ref_col="harmonized_celltype_fine"

##### Load reference and query datasets ####
reference="~/immune-cell-signature-discovery-classification-paper/data/Hao/pbmc_hao_ref_up.qs"
query="./data/Greenland/Greenland_all_sobj_RNA_ADT.qs"

# Load query data
query <- qread(query)

# Load reference data
reference <- qread(reference)

# Group 
# Minor cell type annotations 
reference$harmonized_celltype_fine=reference$celltype.l2
table(reference$harmonized_celltype_fine)

reference$harmonized_celltype_fine=gsub("ASDC","DC",reference$harmonized_celltype_fine)
reference$harmonized_celltype_fine=gsub("cDC1","DC",reference$harmonized_celltype_fine)
reference$harmonized_celltype_fine=gsub("cDC2","DC",reference$harmonized_celltype_fine)


reference$harmonized_celltype_fine=gsub("B intermediate","B",reference$harmonized_celltype_fine)
reference$harmonized_celltype_fine=gsub("B memory","B",reference$harmonized_celltype_fine)
reference$harmonized_celltype_fine=gsub("B naive","B",reference$harmonized_celltype_fine)

reference$harmonized_celltype_fine=gsub("CD14 Mono","Mono",reference$harmonized_celltype_fine)
reference$harmonized_celltype_fine=gsub("CD16 Mono","Mono",reference$harmonized_celltype_fine)

reference$harmonized_celltype_fine=gsub("Plasmablast","Plasma",reference$harmonized_celltype_fine)

reference$harmonized_celltype_fine=gsub("CD4 CTL","TCD4 CTL",reference$harmonized_celltype_fine)
reference$harmonized_celltype_fine=gsub("CD4 Naive","TCD4 naive",reference$harmonized_celltype_fine)

reference$harmonized_celltype_fine=gsub("CD4 TCM","TCD4 TCM",reference$harmonized_celltype_fine)
reference$harmonized_celltype_fine=gsub("CD4 TEM","TCD4 TEM",reference$harmonized_celltype_fine)

reference$harmonized_celltype_fine=gsub("CD8 Naive","TCD8 naive",reference$harmonized_celltype_fine)
reference$harmonized_celltype_fine=gsub("CD8 TCM","TCD8 TCM",reference$harmonized_celltype_fine)
reference$harmonized_celltype_fine=gsub("CD8 TEM","TCD8 TEM",reference$harmonized_celltype_fine)

reference$harmonized_celltype_fine=gsub("dnT","Tunconv dnT",reference$harmonized_celltype_fine)
reference$harmonized_celltype_fine=gsub("gdT","Tunconv gdT",reference$harmonized_celltype_fine)
reference$harmonized_celltype_fine=gsub("MAIT","Tunconv MAIT",reference$harmonized_celltype_fine)

reference$harmonized_celltype_fine=gsub("Treg","TCD4 reg",reference$harmonized_celltype_fine)

reference$harmonized_celltype_fine=gsub("NK_CD56bright","NK",reference$harmonized_celltype_fine)


reference=reference[,!reference$harmonized_celltype_fine %in% c("CD8 Proliferating",
                                                               "CD4 Proliferating",
                                                               "NK Proliferating")]
table(reference$harmonized_celltype_fine)

########## Prepare the gene signatures ############
#### Our genes ####
our_genes= read.table("./published_signatures/immdisc_aybey_final_list.tsv",sep = "\t",
                      stringsAsFactors = F,header = T)
table(our_genes$Annotation)

select=c(                    "Plasma"   ,
                             "NK"                   ,
                             "T CD8"              ,
                             "B"             ,
                             "Monocytes"               ,
                             "DC",
                             "T CD4"   ,
                             "pDC"                           
)

our_genes= our_genes %>%
  subset(Annotation %in%select)
our_genes= our_genes$genes
length(our_genes) #167

print("Geneset ready")

########## RF training and prediction ############

print("RF training and prediction starts")

RF_signatures(geneset="ourgenes",
              signature_genes = our_genes,
              assay = "RNA",
              reference =reference ,
              ref_col=ref_col,
              query =query,
              ref_name=ref_name,
              query_name=query_name,
              save_dir=save_dir)

print("Prediction model and results are saved")
print("RF training and prediction ends")

