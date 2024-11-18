# Preprocess Fujiwara discovery dataset
# Obtained from GEO
# FPKM
save_dir="./data/IFN_bulk/Fujiwara/" 

# Expression data

h_all <- read.table(paste0(save_dir,"/GSE120844_FPKM_table-4.txt"), 
                    header = T,row.names = 1)

dim(h_all)

# Meta data
df=data.frame(Stim=c(rep("untreated",2),
                     rep("IFNg",2),
                     rep("untreated",2),
                     rep("IFNg",2)
),
treatment=c(rep("wt",2),
            rep("wt",2),
            rep("Pom",2),
            rep("Pom",2)
))

# Only take control and IFNg treated samples 

h_all=h_all[,1:4]
df=df[1:4,]
rownames(df)=colnames(h_all)

# Pre-filtering low counts for tpm or fpkm 
keep_genes=rownames(h_all)[rowSums(h_all) > 0]
h_all=h_all[keep_genes,]
dim(h_all)

fuji_list=list(h_all,df)
saveRDS(fuji_list,paste0(save_dir,"fujiwara.RDS"))
