# obtain protein exp for Greenland
# obtain protein markers 
import rpy2.robjects as robjects


# Define your file paths
import os

# Get the absolute path
file_path = os.path.expanduser('~/IFN_DISC/data/Greenland/GSE181897_concat_4_raw.h5ad')

# Read the H5AD file using scanpy
import scanpy as sc
# Check if the file exists
if os.path.exists(file_path):
    adata = sc.read(file_path)
else:
    print(f"File not found: {file_path}")

# Genes and protein markers dataframe
gene_prot_df=adata.var

gene_prot_df['genome'].str.contains("BD99AbSeq").sum() #96 protein 

gene_prot_df.to_csv('Greenland_protein_gene_meta.csv')
