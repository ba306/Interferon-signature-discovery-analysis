# Greenland dataset preprocess

library(Seurat)
library(SeuratDisk)
library(qs)
library(dplyr)

save_dir="data/Greenland/" 

#  GSE181897_concat_4.raw.h5ad from GEO
# Convert h5ad to h5Seurat
Convert(paste0(save_dir,"GSE181897_concat_4.raw.h5ad"),
        paste0(save_dir,"GSE181897_concat_4_raw.h5seurat"))

# Load h5Seurat 
seuratObject <- LoadH5Seurat("data/Greenland/GSE181897_concat_4_raw.h5seurat")
mat=seuratObject@assays$RNA@data
meta=seuratObject@meta.data

# Separate ADT and GEX data
# loading converted the protein names | and _ into -
# gene data obtained Greenland_protein_gene_meta.py
gene_protein_meta_df=read.csv(paste0(save_dir,"Greenland_protein_gene_meta.csv"),row.names = 1)
protein_names=gene_protein_meta_df[grepl("BD99AbSeq",gene_protein_meta_df$genome),] %>%rownames()
protein_names=gsub("|","-",protein_names,fixed = T)
protein_names=gsub("_","-",protein_names,fixed = T)

# Gene matrix GEX
mat_RNA=mat[!rownames(mat)%in% protein_names,]

# Protein matrix ADT
mat_ADT=mat[rownames(mat)%in% protein_names,]

# Create seurat object 
pbmc_f=CreateSeuratObject(counts = mat_RNA,meta.data = meta)


# QC
pbmc_f[["percent.mt"]] <- PercentageFeatureSet(pbmc_f, pattern = "^MT-")
VlnPlot(pbmc_f,features = c("nFeature_RNA"), ncol = 1,
        y.max = 5000)
VlnPlot(pbmc_f,features = c("percent.mt"), ncol = 1,
        y.max = 8)

pbmc_f <- subset(pbmc_f, subset = nFeature_RNA > 200 & 
                   nFeature_RNA < 4000 & percent.mt < 5)

# create a new assay to store ADT information
adt_assay <- CreateAssay5Object(counts = mat_ADT)

pbmc_f[["ADT"]] <- CreateAssayObject(counts = mat_ADT)

# qsave(pbmc_f,paste0(save_dir,
#                "Sun_pbmc_SeuratRobj_RNA_ADT.qs"))
# 
# pbmc_f=qread(paste0(save_dir,
#                       "Sun_pbmc_SeuratRobj_RNA_ADT.RDS"))


# Normalize data 
pbmc_f=NormalizeData(pbmc_f,assay = "RNA",normalization.method = "LogNormalize")
pbmc_f=NormalizeData(pbmc_f,assay = "ADT",normalization.method = "CLR")

# Remove other conditions than IFN and TNFa
# Keep TNFa as negative control 
pbmc_f=subset(pbmc_f,subset = cond!= "0" )
pbmc_f=subset(pbmc_f,subset = cond!= "R" )
pbmc_f=subset(pbmc_f,subset = cond!= "P" )

# Annotate conditions
pbmc_f$cond=gsub("A","TNF-a",pbmc_f$cond)
pbmc_f$cond=gsub("B","IFN-b",pbmc_f$cond)
pbmc_f$cond=gsub("G","IFN-g",pbmc_f$cond)
pbmc_f$cond=gsub("C","Control",pbmc_f$cond)

# Remove not cell type related cells or smaller subpopulations (n<3) based on author annotations
pbmc_f=subset(pbmc_f,subset = ct2!= "M_cDC_PMA/I" )
pbmc_f=subset(pbmc_f,subset = ct2!= "Mitotic" )

# Annotate author cell type annotations
pbmc_f$ct2=gsub("ncM","CD16 Monocyte",pbmc_f$ct2)
pbmc_f$ct2=gsub("cM","CD14 Monocyte",pbmc_f$ct2)
pbmc_f$ct2=gsub("PB","Plasma",pbmc_f$ct2)

# Major cell type annotations 
pbmc_f$cell_types_group=pbmc_f$ct2
pbmc_f$cell_types_group=gsub("B_Mem","B",pbmc_f$cell_types_group)
pbmc_f$cell_types_group=gsub("B_Naive","B",pbmc_f$cell_types_group)
pbmc_f$cell_types_group=gsub("Plasma","B",pbmc_f$cell_types_group)

pbmc_f$cell_types_group=gsub("T4_Mem","T CD4",pbmc_f$cell_types_group)
pbmc_f$cell_types_group=gsub("T4_Naive","T CD4",pbmc_f$cell_types_group)
pbmc_f$cell_types_group=gsub("T_Tox","T_ox",pbmc_f$cell_types_group)
pbmc_f$cell_types_group=gsub("T8_Naive","T CD8",pbmc_f$cell_types_group)

pbmc_f$cell_types_group=gsub("CD14 Monocyte","Myeloid",pbmc_f$cell_types_group)
pbmc_f$cell_types_group=gsub("CD16 Monocyte","Myeloid",pbmc_f$cell_types_group)
pbmc_f$cell_types_group=gsub("cDC","Myeloid",pbmc_f$cell_types_group)
pbmc_f$cell_types_group=gsub("pDC","Myeloid",pbmc_f$cell_types_group)

qsave(pbmc_f,file=paste0(save_dir,"Greenland_all_sobj_RNA_ADT.qs"))

