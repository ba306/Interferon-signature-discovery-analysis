library(readxl)
# Create T CD8 published signatures gene list

save_dir="published_signatures/"

# load T CD8 published signatures 
# activated T cell signatures
published_cell=read.table(file ="published_signatures/all_immune_cell_signatures_curated.txt",sep = "\t",header = T,row.names = 1)

unique(published_cell$Cell_type) %>% sort()

t_cell_names=c(
  "T_CD8"              ,    "T_CD8_activated"     ,   "T_CD8_central_memory" ,  "T_CD8_effector_memory" 
)

published_cell=published_cell[published_cell$Cell_type %in% t_cell_names,]

#published_cell$Gene=unfactor(published_cell$Gene)
published_cell$source_2=paste0(published_cell$Cell_type,".",published_cell$source)
pub_cell=foreach(i=unique(published_cell$source_2))%do% {
  print(i)
  published_cell[published_cell$source_2==i,"Gene"]
}
names(pub_cell)=unique(published_cell$source_2)

# T CD8 cell type signature
refined_df_max_median=read.table("published_signatures/immdisc_aybey_final_list.tsv",
                                 sep = "\t",header = T)
refined_df_max_median=refined_df_max_median[refined_df_max_median$Annotation %in% 
                                              c("T CD8"),]
pub_cell[["T_CD8.Aybey"]]=refined_df_max_median$genes

# Add Nieto T CD8 signature to the signature list
Nieto= read_excel(path = "published_signatures/nieto_immune_2021.xlsx",sheet = 1)

selected=c(     "CD8 effector memory"  ,"CD8 cytotoxic")
nieto_genes=Nieto[Nieto$`Cell type` %in% selected,] %>%.$Markers %>%
  strsplit(", ")

pub_cell[["T_CD8_cytotoxic.Nieto"]] = nieto_genes[[1]]
pub_cell[["T_CD8_eff_memory.Nieto"]] = nieto_genes[[2]]

saveRDS(pub_cell,paste0(save_dir,"T_CD8_published_signatures_list.RDS"))


