library(GEOquery)
library(dplyr)
library(foreach)
library(hgu133plus2.db)
source( "functions/CScores.R")
Sys.setenv(VROOM_CONNECTION_SIZE = 500000)

# GSE121239
# Calculate coherence score, correlation of IFN signatures to 
# SLEDAI for SLE baseline patients

save_dir="SLE_analysis/"

# Download data
gse <- getGEO("GSE121239", GSEMatrix = TRUE)

# Meta data
meta=pData(gse[[1]])
meta$`sledai:ch1` %>%table()
meta$`disease state:ch1` %>%table()

# Select only SLE

meta_1=meta[grepl("v1", meta$title) & meta$`disease state:ch1`=="Systemic Lupus Erythematosus",]

# Expression matrix - rma values in log2 scale
mtx=gse[["GSE121239_series_matrix.txt.gz"]]@assayData[["exprs"]]

mtx=mtx[,rownames(meta_1)]

# Convert into gene symbols
genes_matrix <- gsub("PM_", "", rownames(mtx))

anno <- do.call(cbind, lapply(c( "SYMBOL"),
                              function(x) mapIds(hgu133plus2.db,
          keys=genes_matrix, column = x, keytype = "PROBEID")))

mtx_f=mtx[!is.na(anno[,1]),]
rownames(mtx_f)=anno[!is.na(anno[,1])]

# zscale
mtx_scaled=scale(t(mtx_f))

# all IFN signatures
rosetta_new=readRDS(file ="published_signatures/IFN_signature_list_all_valid.RDS")

meta_plot=meta[rownames(mtx_scaled),]
meta_plot$sledai=meta_plot$`sledai:ch1`
meta_plot=meta_plot[,c("sledai"),drop=F]

save(meta_plot,mtx_scaled,file = paste0(save_dir,"GSE121239_mtx_meta.RData"))
load(file = paste0(save_dir,"GSE121239_mtx_meta.RData"))


#### Mean Signature Scores ####
for(i in seq_along(rosetta_new)){
  genes=intersect(rosetta_new[[i]],colnames(mtx_scaled))
  meta_plot[,names(rosetta_new)[i]]= rowMeans(mtx_scaled[,genes])
}
### Coherence score ####
# add coherence score and number of genes

cor_vector_df=foreach(i=colnames(meta_plot)[-1],.combine=rbind) %do%{
  a=cor.test(meta_plot[,i],as.numeric(meta_plot$sledai), 
             method = "spearman")
  data.frame(signature=i,rho=round(a$estimate[1],2))
}

cs_scores=foreach(i=seq_along(rosetta_new),.combine = rbind) %do% {
  print(i)
  gen= intersect(colnames(mtx_scaled),rosetta_new[[i]])
  pear=cor( mtx_scaled[,gen], method = "pearson")
  score=coherence_score_val(gen, pear)
  round(score,digits = 2)
}

cor_vector_df$signature=cor_vector_df$signature
cor_vector_df$cscore=cs_scores

cor_vector_df$no_genes=lapply(rosetta_new, length) %>%unlist() %>% unname()


cor_vector_df=cor_vector_df %>%arrange(desc(rho))

write.table(cor_vector_df,
            paste0(save_dir,"GSE121239_spearman_cor_ifn_scores_vsledai.txt"),sep="\t", row.names = F)

