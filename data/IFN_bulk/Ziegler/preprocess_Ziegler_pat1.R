# Preprocess Ziegler 
# use only one donor Donor 1 in discovery 
# raw counts received from the authors
# Deseq normalized

library(foreach)
library(dplyr)
library(DESeq2)
source( "./functions/excel_to_genesymbol_ziegler.R") # converting Excel converted months to gene symbols

save_dir="./data/IFN_bulk/Ziegler/" 

# Read expression counts
h1= read.table(paste0(save_dir,"Human1_Basal_Pops_Counts.txt"),
               sep = "\t",stringsAsFactors = F,header = T)

# Take mean of duplicated transcripts
h1_unduplicated=h1 %>%
  group_by(GENE) %>%
  summarise_all(mean)

h1_unduplicated=as.data.frame(h1_unduplicated)
rownames(h1_unduplicated)=h1_unduplicated$GENE

# Remove gene symbol column
h1_unduplicated=h1_unduplicated[,-1]

# detected excel conversion from the original data 
# convert to Gene Symbols!
rownames(h1_unduplicated)=excel_to_genesymbol(h1_unduplicated)

# Meta data
h1_meta= read.table(paste0(save_dir,"Human1 Basal Pops MetaData.txt"),
                    sep = "\t",stringsAsFactors = F,header = T,row.names = 1)

# Sample IDs
h1_meta$rows=rownames(h1_meta)

# Group by stim and within dose-groups

h1_meta=h1_meta %>%
  group_by(Stim) %>%
  arrange(Stim, Dose)

h1_meta=as.data.frame(h1_meta)
rownames(h1_meta)=h1_meta$rows

# Study name
h1_meta$study=rep("ziegler",nrow(h1_meta))

# Take relevant columns
h1_meta=h1_meta[,c("Stim","rows","study","Dose")]

# Remove IL17a and IL4 samples
take_patients=rownames(h1_meta[-which(h1_meta$Stim %in% c("IL17A","IL4")),])

h_all=h1_unduplicated[,take_patients]
h1_meta=h1_meta[take_patients,]

# Donor info 
h1_meta$pat=gsub("COVID_BasalStim_(.*)_.*_.*",h1_meta$rows,replacement = "\\1")

# Expression values as integers 
h_all=as.matrix(h_all)
storage.mode(h_all)="integer"

## Unify IFN annotations

h1_meta$Stim=gsub("IFNA","IFNa",h1_meta$Stim)
h1_meta$Stim=gsub("IFNG","IFNg",h1_meta$Stim)
h1_meta$Stim=gsub("Untreated","untreated",h1_meta$Stim)


# Deseq2
dds <- DESeqDataSetFromMatrix(countData = h_all,
                              colData = h1_meta,
                              design = ~ 1)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- estimateSizeFactors(dds)

dds_counts=counts(dds,normalized=T)

data_list=list(dds_counts,h1_meta)
saveRDS(data_list,paste0(save_dir,"ziegler_donor1.RDS"))


