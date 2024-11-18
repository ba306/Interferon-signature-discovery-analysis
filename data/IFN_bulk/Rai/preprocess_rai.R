# Rai discovery dataset
#https://shiny.ilincs.org/grein
#rai downloaded raw processed from grein

library(dplyr)
library("DESeq2")


save_dir="./data/IFN_bulk/Rai/"


# Expression matrix

rai_raw=read.csv(paste0(save_dir,
                "GSE74863_GeneLevel_Raw_data.csv"),
                header = T)

# Remove Ensembl IDs
rai_raw=rai_raw[,-1]

# Take mean of duplicated transcripts
rai_raw_unduplicated=rai_raw %>%
  group_by(gene_symbol) %>%
  summarise_all(mean)

rai_raw_unduplicated=as.data.frame(rai_raw_unduplicated)
rownames(rai_raw_unduplicated)=rai_raw_unduplicated$gene_symbol
# Remove gene symbol column
data=rai_raw_unduplicated[,-1]


#round to integer as count
data=round(data)

# Meta data
meta=read.csv(paste0(save_dir,
                     "GSE74863_full_metadata.csv"),
              header = T)
rownames(meta)=meta$geo_accession

time=gsub("(.*)_.*_([0-9]*)hrs_Replicate(.*)",meta$title,replacement = "\\2")
treatment=ifelse(grepl("LUC",meta$titl),"wt","shHIRA")
Stim=ifelse(grepl("0hrs",meta$titl),"untreated","IFNb")

meta_df=data.frame(time=time,
                   Stim=Stim,
                   treatment=treatment)
rownames(meta_df)=rownames(meta)

# Take only wild type experiments
df_wt=meta_df[meta_df$treatment=="wt",]
dat_wt=data[,rownames(df_wt)]


# Normalize using DESEQ2
dds <- DESeqDataSetFromMatrix(countData = dat_wt,
                              colData = df_wt,
                              design = ~ 1)
# Remove low expressed genes 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- estimateSizeFactors(dds)

dds_count=counts(dds,normalized=T)

data_list=list(dds_count,df_wt)
saveRDS(data_list,paste0(save_dir,"rai.RDS"))
