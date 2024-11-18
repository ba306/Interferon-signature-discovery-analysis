#Preprocess and save all validation bulk IFN RNA-seq datasets

library(foreach)
library(dplyr)
library("DESeq2")
source( "./functions/excel_to_genesymbol_ziegler.R") # converting Excel converted months to gene symbols


save_dir="./data/IFN_bulk/"

data_dir="./data/IFN_bulk/Validation_datasets/"

# Function to remove low expressed genes
remove_low=function(data){
  rownames(data)[rowSums(data) > 0]
}

###### Ziegler BEAS cell lines ####
# downloaded from GEO GSE148829
# TPM

# Expression matrix
h_beas= read.table(paste0(data_dir,"GSE148829_BEAS_Basal_Pops_TPM.txt"),
                  sep = "\t",stringsAsFactors = F,header = T,row.names = 1)
# remove low expressed genes
h_beas=h_beas[remove_low(h_beas),]

# log2(TPM+1)
h_beas=log2(h_beas+1)

# detected excel conversion from the original data 
# convert to Gene Symbols!
rownames(h_beas)=excel_to_genesymbol(h_beas)

# Meta data
df_beas= read.table(paste0(data_dir,"BEAS Basal Pops MetaData.txt"),
                 sep = "\t",stringsAsFactors = F,header = T,row.names = 1)

df_beas$rows=rownames(df_beas)

# Group by stim and within dose-groups

df_beas=df_beas %>%
  group_by(Stim) %>%
  arrange(Stim, Dose)
df_beas=as.data.frame(df_beas)
rownames(df_beas)=df_beas$rows

# Study name
df_beas$study=rep("ziegler",nrow(df_beas))

# Take relevant columns
df_beas=df_beas[,c("Stim","rows","study","Dose")]

# Remove IL17a and IL4 samples
take_samples_beas=rownames(df_beas[-which(df_beas$Stim %in% c("IL17A","IL4")),])

h_beas=h_beas[,take_samples_beas]
df_beas=df_beas[take_samples_beas,]

# # For plotting purpose add untreated to both stimulation as 0 time point
# 
# add_untreated_a_beas=df_beas[df_beas$Stim=="Untreated",]
# add_untreated_a_beas$Stim=gsub("Untreated","IFNA",add_untreated_a_beas$Stim)
# 
# add_untreated_g_beas=df_beas[df_beas$Stim=="Untreated",]
# add_untreated_g_beas$Stim=gsub("Untreated","IFNG",add_untreated_g_beas$Stim)
# 
# df_beas=df_beas[df_beas$Stim!="Untreated",]
# df_beas=rbind(df_beas,add_untreated_a_beas,add_untreated_g_beas)

ziegler_beas_list=list(h_beas,df_beas)


##### Devlin ######
# downloaded from GEO GSE145647
#normalized by Deseq2 

# Expression matrix
devlin_pbmc_bulk= read.table(paste0(data_dir,"whole_blood_normalized_counts.txt"),
                             sep = "\t",stringsAsFactors = F,header = T,row.names = 1)

# log2(norm. exp+1)
devlin_pbmc_bulk=log2(devlin_pbmc_bulk+1)

# Meta data
# Stimulation
stim=unlist(lapply(strsplit(colnames(devlin_pbmc_bulk),"_"),'[',2))

# remove il4, il3, il10 and il13 treatments
devlin_pbmc_bulk=devlin_pbmc_bulk[,-which(stim %in% c("IL4","IL3","IL10","IL13"))]

stim=unlist(lapply(strsplit(colnames(devlin_pbmc_bulk),"_"),'[',2))

# Donor, replicates
cell_no=unlist(lapply(strsplit(colnames(devlin_pbmc_bulk),"_"),'[',1))

df_devlin=data.frame(Stim=stim,rows=colnames(devlin_pbmc_bulk),
                     study="devlin",cell_no=cell_no)
rownames(df_devlin)=df_devlin$rows

devlin_list=list(devlin_pbmc_bulk,df_devlin)

##### Ziegler donor 2 #########
# downloaded from GEO GSE148829
# TPM

# Expression matrix
h2= read.table(paste0(data_dir,"GSE148829_Human2_Basal_Pops_TPM.txt"),
               sep = "\t",stringsAsFactors = F,header = T,row.names = 1)

# remove low expressed genes
h2=h2[remove_low(h2),]

# log2(TPM+1)
h2=log2(h2+1)

# detected excel conversion from the original data 
# convert to Gene Symbols!
rownames(h2)=excel_to_genesymbol(h2)

# Meta data
h2_meta= read.table(paste0(data_dir,"Human2 Basal Pops MetaData.txt"),
                    sep = "\t",stringsAsFactors = F,header = T,row.names = 1)

h2_meta$rows=rownames(h2_meta)

# Group by stim and within dose-groups

h2_meta=h2_meta %>%
  group_by(Stim) %>%
  arrange(Stim, Dose)
h2_meta=as.data.frame(h2_meta)
rownames(h2_meta)=h2_meta$rows

# Study name

h2_meta$study=rep("ziegler",nrow(h2_meta))

# Take relevant columns

h2_meta=h2_meta[,c("Stim","rows","study","Dose")]

# Remove IL17a and IL4 samples

take_samples_h2=rownames(h2_meta[-which(h2_meta$Stim %in% c("IL17A","IL4")),])

h2=h2[,take_samples_h2]
h2_meta=h2_meta[take_samples_h2,]

# Add donor number
h2_meta$pat=gsub("COVID_BasalStim_(.*)_.*_.*",h2_meta$rows,replacement = "\\1")

# For plotting purpose add untreated to both stimulation as 0 time point

# add_untreated_a_h2=h2_meta[h2_meta$Stim=="Untreated",]
# add_untreated_a_h2$Stim=gsub("Untreated","IFNA",add_untreated_a_h2$Stim)
# 
# add_untreated_g_h2=h2_meta[h2_meta$Stim=="Untreated",]
# add_untreated_g_h2$Stim=gsub("Untreated","IFNG",add_untreated_g_h2$Stim)
# 
# h2_meta=h2_meta[h2_meta$Stim!="Untreated",]
# h2_meta=rbind(h2_meta,add_untreated_a_h2,add_untreated_g_h2)

ziegler_donor2_list=list(h2,h2_meta)


#### Lee 2021 preprocess #######
# Downloaded from GEO GSE161664
# Expression matrices
path_list= list.files(path = paste0(data_dir,"Lee/"),full.names = T,pattern = ".txt")

samples=gsub("(.*)SAEC_(.*).txt",path_list,replacement = "\\2")

# Read expression matrices
lee=foreach(i=seq_along(path_list),.combine = cbind) %do% {
  data= read.table(path_list[i],
                   sep = "\t",stringsAsFactors = F,header = T,row.names = 1)
  colnames(data)=samples[i]
  data
}

# Remove gene names without features, ambiguous reads and not unique aligments 

remove_rows=c("__no_feature" , "__ambiguous" ,    "__alignment_not_unique")
lee=lee[-which(rownames(lee) %in% remove_rows),]

# Meta data
stim=unlist(lapply(strsplit(colnames(lee),"_"),'[',1))

coldata=data.frame(Stim=stim,rows=colnames(lee),study="leeowski")
rownames(coldata)=colnames(lee)

coldata$rep=unlist(lapply(strsplit(coldata$rows,"_"),'[',2))

# Deseq2
dds <- DESeqDataSetFromMatrix(countData = lee,
                              colData = coldata,
                              design = ~ 1)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- estimateSizeFactors(dds)
lee_norm<- counts(dds, normalized=TRUE)

# log2(norm. exp+1)
lee_norm=log2(lee_norm+1)

lee_norm=as.data.frame(lee_norm)
lee_list=list(lee_norm,coldata)

save(ziegler_beas_list,
          devlin_list,
          ziegler_donor2_list,
          lee_list,file = paste0(save_dir,"validation_exp_meta_rdata.RData"))


