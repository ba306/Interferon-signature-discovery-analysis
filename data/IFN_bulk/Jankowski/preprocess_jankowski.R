# Preprocess Jankowski 
# raw counts
# Deseq normalized
# Obtained from GEO

library(foreach)
library(dplyr)
library(DESeq2)

save_dir="./data/IFN_bulk/Jankowski/" 

# Read files
path_list= list.files(path = save_dir,full.names = T,pattern = ".txt")

# Sample names extracted
samples=gsub("(.*)HPPT_(.*).txt",path_list,replacement = "\\2")

# Expression matrix 
jank=foreach(i=seq_along(path_list),.combine = cbind) %do% {
  data= read.table(path_list[i],
                   sep = "\t",stringsAsFactors = F,header = T,row.names = 1)
  colnames(data)=samples[i]
  data
}

# Remove gene names without features, ambiguous reads and not unique aligments 

remove_rows=c("__no_feature" , "__ambiguous" ,    "__alignment_not_unique")
jank=jank[-which(rownames(jank) %in% remove_rows),]


# Meta data

# Stimulation
stim=unlist(lapply(strsplit(colnames(jank),"_"),'[',1))

# Study name
coldata=data.frame(Stim=stim,rows=colnames(jank),study="jankowski")
rownames(coldata)=colnames(jank)

# Replicates numbers
coldata$rep=unlist(lapply(strsplit(coldata$rows,"_"),'[',2))

# Remove IL1b treated samples
jank=jank[,-which(coldata$Stim %in% c("IL1b"))]
coldata=coldata[colnames(jank),]

# Deseq2

dds <- DESeqDataSetFromMatrix(countData = jank,
                              colData = coldata,
                              design = ~ 1)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- estimateSizeFactors(dds)

dds_counts=counts(dds,normalized=T)

data_list=list(dds_counts,coldata)

saveRDS(data_list,paste0(save_dir,"jankowski.RDS"))

