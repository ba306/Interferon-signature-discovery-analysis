library(foreach)
library(dplyr)
library(varhandle)
save_dir="SLE_analysis/"

# Create final ranked list (by correlation) with coherence scores,
# spearman corr scores and number of genes

txt_files=list.files(save_dir,recursive = F,
  pattern = "spearman_cor_ifn_scores_vsledai.txt",
  full.names = T)

gses=gsub("SLE_analysis//(.*)_spearman_cor_ifn_scores_vsledai.txt","\\1",txt_files)

rank_df=foreach(i=seq_along(txt_files),.combine=rbind) %do% {
cor_vector_df=read.table(txt_files[i],header = T,sep = "\t")
cor_vector_df$rank= rank(cor_vector_df$rho)
cor_vector_df
}

# Rank 
sumrank_df=rank_df %>% group_by(signature) %>% summarise(sum(rank)) %>% arrange(desc(`sum(rank)`))
sign_order=sumrank_df$signature

final_tbl=data.frame(matrix(NA,nrow = length(sign_order),ncol =6),
row.names = sign_order)

for(i in seq_along(txt_files)) {
cor_vector_df=read.table(txt_files[i],header = T,sep = "\t",row.names = 1)
cor_vector_df=cor_vector_df[sign_order,]
# rownames(cor_vector_df)=paste0(rownames(cor_vector_df)," (",cor_vector_df$no_genes,")")
rownames(final_tbl)=paste0(rownames(cor_vector_df)," (",cor_vector_df$no_genes,")")
final_tbl[,i]=cor_vector_df$rho
colnames(final_tbl)[i]=paste0(gses[i],"_","rho")
final_tbl[,i+3]=cor_vector_df$cscore
colnames(final_tbl)[i+3]=paste0(gses[i],"_","CS")
}
final_tbl

write.table(final_tbl,
paste0(save_dir,"Final_table_spearman_cor_ifn_scores_vsledai.txt"),
sep="\t", row.names = T)

