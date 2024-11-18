# Colli discovery dataset
# tpm processed downloaded from
#https://www.diabetesepigenome.org/search/?type=Annotation&biosample_term_name=islet+of+Langerhans&annotation_type=gene+expression&annotation_type=variant+allelic+effects&annotation_type_category=RNA-seq&month_released=February%2C+2021

library(dplyr)
library(foreach)
library(data.table)
save_dir="./data/IFN_bulk/Colli/" 

files=list.files(save_dir,full.names = T,pattern = ".tsv")


# Compile TPM expression matrices from each sample 

colli_all_exp=foreach(j=seq_along(files),.combine = cbind)%do%{
print(j)
sample=gsub(pattern = ".*Colli//(.*)_islet.*",x = files[j],
            replacement = "\\1")

colli_tpm=fread(files[j],header = T,sep = "\t")

exp_mat=foreach(i=seq_along(rownames(colli_tpm)),.combine = rbind)%do% {
  print(i)
unlist(strsplit(colli_tpm$TPM_for_all_samples[i],","))
}
colnames(exp_mat)=paste0(sample,"_",1:ncol(exp_mat))
rownames(exp_mat)=colli_tpm$gene
exp_mat
}

storage.mode(colli_all_exp) <- "numeric"

# Pre-filtering low counts
dim(colli_all_exp)
keep_genes=rownames(colli_all_exp)[rowSums(colli_all_exp) > 0]
length(keep_genes)
colli_all_exp=colli_all_exp[keep_genes,]


rownames(colli_all_exp)[duplicated(rownames(colli_all_exp))]


# Take mean of duplicated transcripts

colli_all_exp_df=colli_all_exp
colli_all_exp_df=unlist(colli_all_exp_df)
colli_all_exp_df=as.data.frame(colli_all_exp_df)


# Take mean of duplicated transcripts
colli_all_exp_df$GENE=rownames(colli_all_exp_df)

colli_all_exp_unduplicated=colli_all_exp_df%>%
  group_by(GENE) %>%
  summarise_all(mean)

colli_all_exp_unduplicated=as.data.frame(colli_all_exp_unduplicated)
rownames(colli_all_exp_unduplicated)=colli_all_exp_unduplicated$GENE
colli_all_exp_unduplicated=colli_all_exp_unduplicated[,!grepl("GENE",colnames(colli_all_exp_unduplicated))]

write.table(colli_all_exp_unduplicated,sep="\t",
            file = paste0(save_dir,"colli_TPM_all.txt"))

colli_all_exp_unduplicated=read.table(sep="\t",
            file = paste0(save_dir,"colli_TPM_all.txt"))


# Add meta data

Stim=ifelse(grepl("ctrl",colnames(colli_all_exp_unduplicated)),"untreated","IFNa")
time=gsub(pattern = "[a-z]*([0-9]*)h_.*",x = colnames(colli_all_exp_unduplicated),
     replacement = "\\1")
colli_df=data.frame(Stim=Stim,time=time)

rownames(colli_df)=colnames(colli_all_exp_unduplicated)


colli_list=list(colli_all_exp_unduplicated,colli_df)
saveRDS(colli_list,paste0(save_dir,"colli.RDS"))
       