library(GEOquery)
library(dplyr)
library(foreach)
library(illuminaHumanv4.db)
source( "functions/CScores.R")
Sys.setenv(VROOM_CONNECTION_SIZE = 50000000)
# GSE65391
# Calculate coherence score, correlation of IFN signatures to 
# SLEDAI for SLE baseline patients

save_dir="SLE_analysis/"

# Download data
gse <- getGEO("GSE65391", GSEMatrix = TRUE)

# Meta data
meta=pData(gse[[1]])

# Expression data - batch normalized data
mtx=gse[["GSE65391_series_matrix.txt.gz"]]@assayData[["exprs"]]

# Take only baseline SLE patients

meta_1=meta[meta$`visit:ch1`==1 & meta$`disease state:ch1` =="SLE",]

# log2 normalization
mtx=log2(mtx+1)
mtx=mtx[,rownames(meta_1)]

# ILM ids to Gene IDs

df=data.frame(Gene=unlist(mget(x = rownames(mtx),envir = illuminaHumanv4SYMBOL)))
df=na.omit(df)
mtx=mtx[rownames(df),]
rownames(mtx)=df$Gene

# zscale 
mtx_scaled=scale(t(mtx))

# all IFN signatures
rosetta_new=readRDS(file ="published_signatures/IFN_signature_list_all_valid.RDS")

meta_plot=meta[rownames(mtx_scaled),]
meta_plot$sledai=meta_plot$`sledai:ch1`
meta_plot$sledai=ifelse(is.na(meta_plot$sledai),0,meta_plot$sledai)
meta_plot=meta_plot[,c("sledai"),drop=F]

save(list = c("meta_plot","mtx_scaled"),file = paste0(save_dir,"GSE65391_mtx_meta.RData"))
load(file = paste0(save_dir,"GSE65391_mtx_meta.RData"))


#### Mean Signature Scores ####
for(i in seq_along(rosetta_new)){
  genes=intersect(rosetta_new[[i]],colnames(mtx_scaled))
  meta_plot[,names(rosetta_new)[i]]= rowMeans(mtx_scaled[,genes])
}

cor_vector_df=foreach(i=colnames(meta_plot)[-1],.combine=rbind) %do%{
  a=cor.test(meta_plot[,i],as.numeric(meta_plot$sledai), 
             method = "spearman")
  data.frame(signature=i,rho=round(a$estimate[1],2))
}

### Coherence score ####
# add coherence score and number of genes

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
            paste0(save_dir,"GSE65391_spearman_cor_ifn_scores_vsledai.txt"),sep="\t", row.names = F)
