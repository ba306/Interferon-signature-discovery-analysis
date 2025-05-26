# ROC curves
# Prediction of IFN signatures in ICI compared with IFN-II Aybey

library(pROC)
library(ProfilerAPI2)
library(dplyr)
library(foreach)
library(ggplot2)
library(rstatix)
library(ggpubr)
library(reshape2)
library(gridExtra)
library(grid)
library(readxl)

# Output directory
data_dir="data/Tumor_bulk/"

# Our IFN signatures
rosetta_new=readRDS(file = "published_signatures/IFN_signature_list_all_valid.RDS")

auc_results = data.frame(Signature=names(rosetta_new)) # Initialize


# All genes to extract
all_genes=unique(unlist(rosetta_new))

#### Datasets and analysis #### 

#### Imvigor210 Bladder ####
# Get expression data + metadata 

# get connection to profiler
api <- ProfilerAPI2::profiler_api(profile = "default")
# Get expression data from Profiler
imvigor <- 
  dplyr::tbl(api$conn,"genentech_imvigor210_rna_wts_gsnap_21q1_prod_jami_0") %>% 
  filter(gene_name %in% all_genes) %>%
  select(subject_id,gene_name,best_confirmed_overall_response,binaryresponse,gene_tpm,
         received_platinum,sample_collected_pre_platinum) %>%
  collect()

#imvigor=saveRDS(paste0(data_dir,"IMVIGOR210_exp_meta_IFNgenes_All.RDS"))
imvigor=readRDS(paste0(data_dir,"IMVIGOR210_exp_meta_IFNgenes_All.RDS"))

# remove samples without response info
imvigor=imvigor[!is.na(imvigor$binaryresponse),]

# Responders: SD/PD; non-responders: CR/PR 
imvigor$binaryresponse=
  ifelse(imvigor$binaryresponse=="SD/PD",
         "NR","R")

# Cohort of platinum treated patients are chosen 
# Save individual plots as lists and compile together at the end

# Initialize AUC results list
imvigor_auc = list()

for (g in seq_along(rosetta_new)) {
  
  print(g)
  
  # Calculate signature scores (as you already did)
  score_Df = imvigor %>%
    filter(received_platinum == "Y") %>%
    filter(gene_name %in% unlist(rosetta_new[g])) %>%
    mutate(gene_logtpm = as.vector(scale(log2(1 + as.numeric(gene_tpm))))) %>%
    group_by(subject_id, binaryresponse) %>%
    dplyr::summarize(score = mean(gene_logtpm, na.rm = TRUE), .groups = "drop")
  
  # ROC analysis
  roc_obj = roc(score_Df$binaryresponse, score_Df$score, levels = c("NR", "R"))
  auc_val = as.numeric(auc(roc_obj))
  
  # Store
  imvigor_auc[[names(rosetta_new)[g]]] = auc_val
}

# Other ICB datasets are
# downloaded from https://github.com/xmuyulab/ims_gene_signature

#### Kim gastric #####

# Expression data
kim_mtx=read.table(paste0(data_dir,"gas_korean_exp_data.csv"),
                   sep = ",",header = T,row.names = 1,check.names = F)
# Clinical data
kim_meta=read.table(paste0(data_dir,"gas_korean_cli_data.csv"),
                    sep = ",",header = T,row.names = 1)
# Binary response 
kim_meta$response=case_when(
  kim_meta$response=="-1" ~ "NR",
  kim_meta$response=="0" ~ "NR",
  kim_meta$response=="1" ~ "R",
  T~NA
)

# Select baseline samples
kim_meta=kim_meta%>%
  dplyr::filter(treatment=="pre")
# Select samples with exp and clinical data
kim_inter_pat=intersect(rownames(kim_meta),colnames(kim_mtx))
kim_meta=kim_meta[kim_inter_pat,]

# select baseline samples in expression data 
# z-scale the expression data
kim_mtx_scaled=scale(t(kim_mtx[,kim_inter_pat]))

# Calculate mean signature scores (signatures) or gene expression (individual genes)
# Add to the meta data

for(j in seq_along(rosetta_new)){
  genes=intersect(rosetta_new[[j]],colnames(kim_mtx_scaled))
  kim_meta[,names(rosetta_new)[j]]= rowMeans(kim_mtx_scaled[,genes,drop=F],na.rm = T)
}

# Generate boxplots
kim_meta_plot_melt=melt(kim_meta)

# AUC 

kim_auc = list()

for (g in seq_along(rosetta_new)) {
  
  print(g)
  
  # Calculate signature scores
  genes = intersect(rosetta_new[[g]], colnames(kim_mtx_scaled))
  
  if (length(genes) > 0) {
    
    score_Df = data.frame(
      subject_id = rownames(kim_meta),
      response = kim_meta$response,
      score = rowMeans(kim_mtx_scaled[, genes, drop = FALSE], na.rm = TRUE)
    )
    
    # ROC analysis
    roc_obj = roc(score_Df$response, score_Df$score, levels = c("NR", "R"), direction = "<")
    auc_val = as.numeric(auc(roc_obj))
    
    # Store
    kim_auc[[names(rosetta_new)[g]]] = auc_val
  } else {
    kim_auc[[names(rosetta_new)[g]]] = NA
  }
}


#### Van Allen melanoma ####
# Expression data
allen_mtx=read.table(paste0(data_dir,"mel_van_exp_data.csv"),
                     sep = ",",header = T,row.names = 1,check.names = F)
# Clinical data
allen_meta=read.table(paste0(data_dir,"mel_van_cli_data.csv"),
                      sep = ",",header = T,row.names = 1)
# Binary response 
allen_meta$response=case_when(
  allen_meta$response=="-1" ~ "NR",
  allen_meta$response=="0" ~ "NR",
  allen_meta$response=="1" ~ "R",
  T~NA
)

# Select baseline samples
allen_meta=allen_meta%>%
  dplyr::filter(treatment=="pre")
# Select samples with exp and clinical data
allen_inter_pat=intersect(rownames(allen_meta),colnames(allen_mtx))
allen_meta=allen_meta[allen_inter_pat,]

# select baseline samples in expression data 
# z-scale the expression data
allen_mtx_scaled=scale(t(allen_mtx[,allen_inter_pat]))

# Calculate mean signature scores (signatures) or gene expression (individual genes)
# Add to the meta data

for(j in seq_along(rosetta_new)){
  genes=intersect(rosetta_new[[j]],colnames(allen_mtx_scaled))
  allen_meta[,names(rosetta_new)[j]]= rowMeans(allen_mtx_scaled[,genes,drop=F],na.rm = T)
}

# AUC

vanallen_auc = list()

# Van Allen dataset
for (g in seq_along(rosetta_new)) {
  
  print(g)
  
  # Calculate signature scores
  genes = intersect(rosetta_new[[g]], colnames(allen_mtx_scaled))
  
  if (length(genes) > 0) {
    
    score_Df = data.frame(
      subject_id = rownames(allen_meta),
      response = allen_meta$response,
      score = rowMeans(allen_mtx_scaled[, genes, drop = FALSE], na.rm = TRUE)
    )
    
    # ROC analysis
    roc_obj = roc(score_Df$response, score_Df$score, levels = c("NR", "R"), direction = "<")
    auc_val = as.numeric(auc(roc_obj))
    
    # Store
    vanallen_auc[[names(rosetta_new)[g]]] = auc_val
  } else {
    vanallen_auc[[names(rosetta_new)[g]]] = NA
  }
}

### Save everything into one dataframe
auc_results = data.frame(
  Signature = names(rosetta_new),
  IMvigor210 = unlist(imvigor_auc),
  Kim = unlist(kim_auc),
  VanAllen = unlist(vanallen_auc)
)

# Calculate the mean AUC across datasets
auc_results$Mean_AUC = rowMeans(auc_results[, c("IMvigor210", "Kim", "VanAllen")], na.rm = TRUE)

# Sort by Mean_AUC in descending order
auc_results = auc_results[order(-auc_results$Mean_AUC), ]

# round
auc_results[,-1]=round(auc_results[,-1],2)

# Save sorted AUC results
write.csv(auc_results, paste0(save_dir, "ICB_IFN_signature_AUCs_ranked.csv"), row.names = FALSE)









