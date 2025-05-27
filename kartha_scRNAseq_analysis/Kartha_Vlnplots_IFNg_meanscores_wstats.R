# Analysis on Kartha dataset IFNg 1h, 6h
# Violin plots for different IFN signature scores
# UMAP based on our cell type gene set followed by showing our IFN signature scores on UMAP

library(Seurat)
library(uwot)
library(ggrepel)
library(dplyr)
library(reshape2)
library(varhandle)
library(qs)
library(ggpubr)
options(ggrepel.max.overlaps = Inf)
set.seed(42)

save_dir="kartha_scRNAseq_analysis/" 

# Data taken from our Front Imm Paper 2023 
# https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2023.1194745/full
data_dir="data/Kartha/" 

### Load IFN gene signatures
rosetta_new=readRDS(file ="published_signatures/IFN_signature_list_all_valid.RDS")

# Only use our IFN-g  signature 
selected_signatures=rosetta_new$IFN_II_Aybey

# Load Greenland data
pbmc=qread(file=paste0(data_dir,"Kartha_sobj.qs"))


# cell type assignments from our previous paper 
cell_annot=read.csv(paste0(data_dir,"rf_Hao_fine_Kartha_ourgenes_164.csv"),row.names = 1)

# group cell type annotation
cell_annot$Predictions=gsub("TCD8 TEM","TCD8 memory",cell_annot$Predictions)
cell_annot$Predictions=gsub("TCD8 TCM","TCD8 memory",cell_annot$Predictions)
cell_annot$Predictions=gsub("TCD4 TEM","TCD4 memory",cell_annot$Predictions)
cell_annot$Predictions=gsub("TCD4 TCM","TCD4 memory",cell_annot$Predictions)

pbmc[['fine_RF_annot']]=cell_annot[colnames(pbmc),"Predictions"]

# Only focus on IFN and control for the demonstration of downstream analysis bias

take=pbmc$StimType=="Control" | pbmc$StimType=="IFN"
pbmc=pbmc[,take]

# Remove ControlGolgiPlug_6h and IFNGolgiPlug_6h
# remove donor 4 since this only has IFNg treatment
remove=pbmc$Condition=="ControlGolgiPlug_6h" | pbmc$Condition=="IFNGolgiPlug_6h" | pbmc$Donor =="Donor4"
pbmc=pbmc[,!remove]

# remove small number of cells ILC and T unconv dnT
# each stimulation less than 10
table(pbmc$fine_RF_annot,pbmc$Condition)

pbmc=pbmc[,pbmc$fine_RF_annot!="Tunconv gdT"]
pbmc=pbmc[,pbmc$fine_RF_annot!="Plasma"]
pbmc=pbmc[,pbmc$fine_RF_annot!="HSPC"]
pbmc=pbmc[,pbmc$fine_RF_annot!="pDC"]
pbmc=pbmc[,pbmc$fine_RF_annot!="Tunconv MAIT"]

# only focus on  IFN genes for the downstream analysis
all_genes=unlist(selected_signatures,use.names = F) %>% unique()
pbmc_IFN=pbmc[all_genes,]

#### IFN scores ####

# scale gene expression
pbmc_IFN=ScaleData(pbmc_IFN,features = rownames(pbmc_IFN))

# Calculate mean signature scores

genes_select=intersect(rownames(pbmc_IFN),selected_signatures)
pbmc_IFN[['IFN_II_Aybey']]=colMeans(pbmc_IFN@assays$RNA@scale.data[genes_select,])
  
# Create dataframe for plotting
df=data.frame(
  stimulation=pbmc_IFN$Condition,
  cell_type_fine=pbmc_IFN$fine_RF_annot,
  
  IFN_II_Aybey=pbmc_IFN@meta.data[,"IFN_II_Aybey"]
)

df_melt=melt(df)

# Predefined parameters
text_size    = 14
signif_size  = 6
bracket_size = 1

# Specify treatment IFNg
df$stimulation=gsub("IFN","IFNg",df$stimulation)

# Prepare data by creating helper columns and converting groups to factors
df2 <- df %>% 
  mutate(
    Time = ifelse(grepl("1h", stimulation), "1h", "6h"),
    Group = ifelse(grepl("Control", stimulation), "Control", "IFNg")
  ) %>%
  mutate(
    Group = factor(Group, levels = c("Control", "IFNg")),
    cell_type_fine = factor(cell_type_fine)
  )

# Create separate datasets for 1h and 6h comparisons
df_1h <- df2 %>% filter(Time == "1h")
df_6h <- df2 %>% filter(Time == "6h")

# Function to generate violin plot with custom t-test
plot_violin_with_stats <- function(data, time_label, ref_group) {
  # For each cell type, perform t-test comparing groups
  pwc <- data %>% 
    group_by(cell_type_fine) %>% 
    t_test(IFN_II_Aybey  ~ Group, 
           ref.group = ref_group)
  
  # If p.adj is missing, use p instead.
  if(!"p.adj" %in% colnames(pwc)){
    pwc <- pwc %>% mutate(p.adj = p)
  }
  
  # Create significance labels based on adjusted p-values
  pwc <- pwc %>% 
    mutate(p.adj.signif = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE ~ "ns"
    ))
  
  pwc <- pwc %>% add_xy_position(x = "cell_type_fine")
  
  # Create violin plot with custom p-value annotations
  p <- ggplot(data, aes(x = cell_type_fine, y = IFN_II_Aybey , fill = Group)) +
    geom_violin(alpha = 0.5, scale = "width", position = position_dodge(width = 1)) +
    ggbeeswarm::geom_quasirandom(shape = 21, size = 2, dodge.width = 1,
                                 color = "black", alpha = 0.5, show.legend = FALSE) +
    ggtitle(time_label) +
    theme_bw() +
    theme(legend.position = "top",
          legend.text = element_text(size = text_size),
          axis.text.x = element_text(size = text_size, angle = 30, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = text_size),
          strip.text.x = element_text(size = text_size)) +
    ylab("Mean signature score") + xlab("Cell types")
  
  # Use inherit.aes = FALSE to avoid inheriting unnecessary mappings
  p <- p + stat_pvalue_manual(
    data = pwc,
    mapping = aes(x = cell_type_fine, y.position = y.position, label = p.adj.signif),
    tip.length = 0.01, hide.ns = F, 
    size = signif_size, bracket.size = bracket_size,
    inherit.aes = FALSE
  )
  
  return(p)
}

# Create the violin plots for 1h and 6h comparisons:
plot_1h <- plot_violin_with_stats(df_1h, "1h", ref_group = "IFNg")
plot_6h <- plot_violin_with_stats(df_6h, "6h", ref_group = "IFNg")

# Display plots together

vln_all_common_3=plot_grid(plot_1h,plot_6h,nrow = 2, align = "h")

pdf(paste0(save_dir,"Kartha_Vlnplots_IFNg_meanscores_wstats.pdf"),
    width = 8,
    height = 8)
vln_all_common_3
dev.off()

