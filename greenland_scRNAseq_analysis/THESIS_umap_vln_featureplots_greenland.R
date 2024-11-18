# Analysis on Greenland dataset
# UMAP based on our cell type gene set followed by showing our IFN signature scores on UMAP
# Violin plots for different IFN signature scores

library(Seurat)
library(uwot)
library(ggrepel)
library(dplyr)
library(reshape2)
library(varhandle)
library(qs)
library(ggpubr)
options(ggrepel.max.overlaps = Inf)
library(emmeans)
library(rstatix)
library(patchwork)

set.seed(42)

save_dir="greenland_scRNAseq_analysis/" 
data_dir="data/Greenland/" 

# Color palette
chr=c("darkblue","brown4","darkgreen","darkorange","darkcyan",
      "darksalmon","deeppink","red2", "goldenrod1", "mediumorchid",
      "#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0","black")

### Load IFN gene signatures
rosetta_new=readRDS(file ="published_signatures/IFN_signature_list_all_valid.RDS")

# Only use IFN-g published signatures and our signatures 
selected_signatures=rosetta_new[c("IFN_I_Aybey","IFN_II_Aybey",
                                  names(rosetta_new)[grepl("ifng",
                                                           names(rosetta_new),ignore.case = T)])]

# Load Greenland data
pbmc=qread(file=paste0(data_dir,"Greenland_all_sobj_RNA_ADT.qs"))



# cell type annototation from Random forest classifier
greenland_annot=read.csv("data/Greenland/rf_Hao_fine_Greenland_ourgenes_165.csv",row.names = 1)
# group cell type annotation
greenland_annot$Predictions=gsub("TCD8 TEM","TCD8 memory",greenland_annot$Predictions)
greenland_annot$Predictions=gsub("TCD8 TCM","TCD8 memory",greenland_annot$Predictions)
greenland_annot$Predictions=gsub("TCD4 TEM","TCD4 memory",greenland_annot$Predictions)
greenland_annot$Predictions=gsub("TCD4 TCM","TCD4 memory",greenland_annot$Predictions)


pbmc[['fine_RF_annot']]=greenland_annot[colnames(pbmc),"Predictions"]

table(pbmc$fine_RF_annot)
table(pbmc$fine_RF_annot,pbmc$cond)
table(pbmc$cell_types_group)

# remove small number of cells ILC and T unconv dnT
# each stimulation less than 10

pbmc=pbmc[,pbmc$fine_RF_annot!="ILC" & 
            pbmc$fine_RF_annot!="Tunconv dnT"]

pbmc=pbmc[,pbmc$fine_RF_annot!="Plasma"]


# Remove - from conditions 
pbmc$cond=gsub("-","",pbmc$cond)

# only focus on  IFN genes for the downstream analysis
all_genes=unlist(rosetta_new,use.names = F) %>% unique()
pbmc_IFN=pbmc[all_genes,] #308 genes and 89498 samples 

#### IFN scores ####

# scale gene expression
pbmc_IFN=ScaleData(pbmc_IFN,features = rownames(pbmc_IFN))

# Calculate mean signature scores
for(i in seq_along(selected_signatures)){
  genes_select=intersect(rownames(pbmc_IFN),selected_signatures[[i]])
  pbmc_IFN[[names(selected_signatures)[i]]]=colMeans(pbmc_IFN@assays$RNA@scale.data[genes_select,])
  
}  



# Create dataframe for plotting
df=data.frame(
  stimulation=pbmc_IFN$cond,
  cell_type_fine=pbmc_IFN$fine_RF_annot,
  
  pbmc_IFN@meta.data[,names(selected_signatures)]
)

df_melt=melt(df)


# pairwise comparisons


text_size=40
signif_size=12
bracket_size=1.5

df_melt$stimulation=gsub("-","",df_melt$stimulation)

unq_variable=unique(df_melt$variable)

vlnmeanscore_stat=foreach(i= unq_variable) %do% {
  print(i)
  if(grepl("IFN_I_",i)==T){
    pwc <-df_melt %>%
      filter(variable==i)%>% 
      group_by(cell_type_fine) %>%
      t_test(value ~ stimulation, p.adjust.method = "bonferroni",ref.group = "IFNb",alternative = "greater")
  } else{
    pwc <-df_melt %>%
      filter(variable==i)%>% 
      group_by(cell_type_fine) %>%
      t_test(value ~ stimulation, p.adjust.method = "bonferroni",ref.group = "IFNg",alternative = "greater")
  }
  pwc
  
  
  # Visualization: box plots with p-values
  pwc <- pwc %>% add_xy_position(x = "cell_type_fine")
  pwc
  
  print(df_melt %>% 
          filter(variable==i)%>%
          ggplot(aes(x=cell_type_fine, y=value
          )) +
          geom_violin(aes(fill = stimulation), trim = FALSE,scale = "width") +
          theme_bw()+  
          scale_fill_brewer(palette="Dark2")+
          facet_wrap(facets = "variable",nrow = 5,scales = "free")+
          theme(legend.position='none',
                legend.text = element_text(size = text_size),
                axis.text.x = element_text(size = text_size,angle = 30, vjust = 1, hjust=1),
                axis.text.y = element_text(size = text_size),
                strip.text.x = element_text(size = text_size))+
          ylab("")+xlab("")+
          stat_pvalue_manual(pwc,label = "p.adj.signif", 
                             tip.length = 0.01,hide.ns = T,size=signif_size,
                             bracket.size = bracket_size)
  )
  
}


# Get legend 
legend <- get_legend(
  df_melt %>% 
    rename(Stimulation = stimulation) %>%
    filter(variable==i)%>%
    ggplot(aes(x=cell_type_fine, y=value
    )) +
    geom_violin(aes(fill = Stimulation), trim = FALSE,scale = "width") +
    theme_bw()+  
    scale_fill_brewer(palette="Dark2")+
    facet_wrap(facets = "variable",nrow = 5,scales = "free")+
    theme(legend.position='top',text = element_text(size = text_size),
          legend.text = element_text(size = text_size))+
    ylab("")+xlab("")+
    stat_pvalue_manual(pwc,label = "p.adj.signif",label.size = 6)) 

plot(legend)

# Add a common y-axis label
common_label <- ggdraw() +
  draw_label("Mean signature score", size = 60, angle = 90,fontface = "plain")

vln_all=wrap_plots(vlnmeanscore_stat,ncol = 2)

# Combine the grid and common label
vln_all_common=plot_grid(common_label,vln_all ,ncol = 2, align = "v",rel_widths = c(0.04,1))
vln_all_common_2=plot_grid(legend,vln_all_common ,nrow = 2, align = "h",rel_heights  = c(0.02,1))


common_label_y <- ggdraw() +
  draw_label("Cell types", size = 60, fontface = "plain")
vln_all_common_3=plot_grid(vln_all_common_2,common_label_y,nrow = 2, align = "h",rel_heights = c(1,0.02))


pdf(paste0(save_dir,"Vlnplots_IFN_meanscores_wstats.pdf"),
    width = 30,
    height = 45)
vln_all_common_3
dev.off()

#### UMAP based on our cell type genes and featureplots for our IFN signatures ####

# Our cell type genes
db_genes= read.table("published_signatures/immdisc_aybey_final_list.tsv",sep = "\t",
                     stringsAsFactors = F,header = T)

select=c(                    "Plasma"   ,
                             "NK"                   ,
                             "T CD8"              ,
                             "B"             ,
                             "Monocytes"               ,
                             "DC",
                             "T CD4"   ,
                             "pDC"                           
)

db_genes_filter= db_genes %>%
  subset(Annotation %in%select)
cellgenes= db_genes_filter$genes

pbmc_f_celltype=pbmc[cellgenes,] #165 genes
dim(pbmc_f_celltype)

pbmc_f_celltype=ScaleData(pbmc_f_celltype,features = rownames(pbmc_f_celltype))

# UMAP based on our cell type genes
pbmc_f_celltype <- RunPCA(pbmc_f_celltype, 
                          features = rownames(pbmc_f_celltype))
ElbowPlot(pbmc_f_celltype)
pbmc_f_celltype <- RunUMAP(pbmc_f_celltype, dims = 1:5)

# um_ourgenes_cell_types showing cell types 
um_ourgenes_cell_types=DimPlot(pbmc, reduction = "umap",group.by="fine_RF_annot",label=T,pt.size = 2,
                               repel=T,cols=chr,label.box = T,label.color = "white",label.size = 14,raster=T)+ 
  theme(legend.position = c(0.8, 0.8),
        text = element_text(size = 40),
        axis.text.x = element_text(size = 40),  # Increase x-axis labels size
        axis.text.y = element_text(size = 40),
        plot.title = element_text(size = 40, face = "bold", hjust = 0.5))+ggtitle("RF cell types")

# um_ourgenes_stim showing treatments

um_ourgenes_stim=DimPlot(pbmc, reduction = "umap",group.by="cond",label=F,pt.size = 2,raster = T,
                         label.size = 18,repel=T)+ 
  theme(legend.position = c(0.02, 0.6),
        text = element_text(size = 40),
        axis.text.x = element_text(size = 40),  # Increase x-axis labels size
        axis.text.y = element_text(size = 40),
        plot.title = element_text(size = 40, face = "bold", hjust = 0.5))+ggtitle("Condition")

fig_up=patchwork::wrap_elements(um_ourgenes_cell_types) +
  patchwork::wrap_elements(um_ourgenes_stim)

# Feature plot for IFN-aybey scores
pbmc_IFN@reductions=reds
fig_down=FeaturePlot(pbmc_IFN, 
                     features = names(selected_signatures)[grepl("aybey", names(selected_signatures), ignore.case = TRUE)], 
                     split.by = "cond", 
                     cols = c("gray", "red"),raster = T,pt.size = 2) & 
  theme_bw() + 
  theme(
    text = element_text(size = 40),
    legend.position = "top",              # Move legend to the top
    legend.direction = "horizontal"         # Arrange the legend vertically
  )

pdf(paste0(save_dir,"UMAP_Featureplot_IFN_aybey.pdf"),width = 40,
    height = 35)

plot_grid(fig_up,NULL,fig_down,nrow=3,rel_heights = c(1,0.1,1))

dev.off()

#### UMAP on IFN genes ####

# How well are different IFN genes are separated?
# Include only control and IFN treatments

pbmc_IFN_woTNF=subset(pbmc_IFN,subset = cond!= "TNF-a" )
pbmc_IFN_woTNF=ScaleData(pbmc_IFN_woTNF,features = rownames(pbmc_IFN_woTNF))

# Categorize IFN signatures --> our IFN genes  vs IFN-I or IFN-II published genes

rosetta_new_df=melt(rosetta_new)

rosetta_new_df$IFN_Type=
  ifelse(grepl("aybey",rosetta_new_df$L1,ignore.case = T),rosetta_new_df$L1,
         ifelse(grepl("IFNg",rosetta_new_df$L1,ignore.case = T),
                "IFN-II","IFN-I"))

# IFNI published genes
IFNI_genes=unique(rosetta_new_df[rosetta_new_df$IFN_Type=="IFN-I","value"])
length(IFNI_genes) #186

# IFNII published genes
IFNII_genes=unique(rosetta_new_df[rosetta_new_df$IFN_Type=="IFN-II","value"])
length(IFNII_genes) #229

# Common genes for IFN-I and IFN-II published genes
common=intersect(IFNI_genes,IFNII_genes)
length(common) #95

# Scaled expression matrix
data_j=as.matrix(pbmc_IFN_woTNF@assays$RNA@scale.data)

# UMAP based on our IFN genes

data_j_i_ii=data_j[rownames(data_j) %in% c(rosetta_new$IFN_I_Aybey,
                                           rosetta_new$IFN_II_Aybey),]

umap_clustergenes= umap(data_j_i_ii,pca=2,verbose = T)

um_clustergenes =as.data.frame(umap_clustergenes)

um_clustergenes$IFNsignature=ifelse(rownames(um_clustergenes)%in%rosetta_new$IFN_I_Aybey, 
                                    "IFN-I Aybey",
                                    ifelse(rownames(um_clustergenes)%in%rosetta_new$IFN_II_Aybey,
                                           "IFN-II Aybey","")
) 

um_clustergenes$gene=rownames(um_clustergenes)
our_umap=ggplot(um_clustergenes,aes(V1,V2,colour=IFNsignature)) +
  geom_point(size=5) + theme_bw()+
  theme( legend.position="top", 
         legend.spacing.x = unit(0.1, 'cm'),
         text = element_text(size = 30),
         axis.text.x = element_text(size = 30),  # Increase x-axis labels size
         axis.text.y = element_text(size = 30))+
  geom_text_repel(aes(label =  as.character (gene)),
                  nudge_x = 0.1, direction = "y",size =10,show.legend = F,force = 10)+
  xlab( "UMAP_1")+ylab("UMAP_2")+ 
  scale_color_discrete(name = "IFN signature")


## umap based on ifn i and ifn ii published gene sets
data_j_i_ii_pub=data_j[rownames(data_j) %in% c(IFNI_genes,IFNII_genes),]

umap_clustergenes_pub= umap(data_j_i_ii_pub,pca=2,verbose = T)

umap_clustergenes_pub =as.data.frame(umap_clustergenes_pub)

umap_clustergenes_pub$IFNsignature=ifelse(rownames(umap_clustergenes_pub)%in%common, 
                                          "IFN-I/II common",
                                          ifelse(rownames(umap_clustergenes_pub)%in%IFNII_genes,
                                                 "IFN-II","IFN-I")
) 

umap_clustergenes_pub$gene=rownames(umap_clustergenes_pub)
published_umap=ggplot(umap_clustergenes_pub,aes(V1,V2,color=IFNsignature)) +
  geom_point(size=5) +  theme_bw()+
  theme(legend.position="top", 
        legend.spacing.x = unit(0.1, 'cm'),
        text = element_text(size = 30),
         axis.text.x = element_text(size = 30),  # Increase x-axis labels size
         axis.text.y = element_text(size = 30))+
  xlab( "UMAP_1")+ylab("UMAP_2")+ 
  scale_color_discrete(name = "IFN signature")

# All together
pdf(paste0(save_dir,"UMAP_Greenland_published_vs_ours.pdf"),
    width=25,
    height=13)
ggarrange(our_umap,published_umap,labels="AUTO",font.label = list(size=40))

dev.off()



