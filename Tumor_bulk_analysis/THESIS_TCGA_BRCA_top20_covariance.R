library(xopdata) # internal package
library(dplyr)
library(foreach)
library(reshape2)
library(readxl)
library("ComplexHeatmap")
library(gridtext)
library(ProfilerAPI2) # internal package
source("functions/RosettaSX_heatmap.R")

# Comparing our IFN signatures against other signatures in a similar analysis done by 
# Kreis et al. in RosettaSX analysis on TCGA BRCA
# Top 20 

# Example cohort TCGA BRCA

# Signatures are extracted from Julian Kreis' analysis- pre-calculated for https://www.rosettasx.com/
# New signatures for T CD8 and IFN-II signatures are added
# T CD8 signature are added from previously curated list incl. our signature
# IFN-II signatures are also added from curated IFN signature list

save_dir="Tumor_bulk_analysis/"

src = xopdata::get_src(pattern="tcga") # connect to the server

#### Load additional signatures ####

# Load IFN-II-IFN-I signatures
rosetta_new=readRDS(file = "published_signatures/IFN_signature_list_all_valid.RDS")

# Take IFN signatures not present in RosettaSX
taken_signatures=names(rosetta_new)[grepl("^IFNg",names(rosetta_new)) &
                                      !grepl("IFNg_Hallmark",names(rosetta_new)) &
                                      !grepl("IFNg_Dummer",names(rosetta_new)) ]

data_list=rosetta_new[c("IFN_I_Aybey",
                        "IFN_II_Aybey",
                        taken_signatures )]

# load T CD8 published signatures 
T_cd8_sign=readRDS(file = "published_signatures/T_CD8_published_signatures_list.RDS")

# Add T CD8 signatures to the IFN-II signature list
data_list=c(data_list,T_cd8_sign)

# Query over analyzed genes only!
all_genes=unname(unlist(data_list)) #411

#### Get data ####
crc_tmm_gexp = merck_tcga_rna_wtstmm_gexp() %>%
  dplyr::filter(id_gene_symbol %in% all_genes,
                id_tcga_cohort %in% c("BRCA")) %>%
  dplyr::select(id_tcga_sample_barcode, id_gene_symbol,
                quant_tmm_tpm) %>%
  dplyr::collect()

crc_tmm_gexp <- saveRDS(crc_tmm_gexp,"data/Tumor_bulk/crc_tmm_gexp.RDS")
crc_tmm_gexp <- readRDS("data/Tumor_bulk/crc_tmm_gexp.RDS")

#### Coherence score ####
# also calculate fraction of genes found/ whole signature
CS_scores_df=foreach::foreach(i=names(data_list),.combine = rbind) %do% {
  print(i)
  
  CS=crc_tmm_gexp%>%
    dplyr::filter(id_gene_symbol %in% data_list[[i]]) %>% 
    group_by(id_gene_symbol) %>%
    dplyr::mutate(quant_tmm_tpm = as.vector(scale(log2(quant_tmm_tpm + 1)))) %>%
    dplyr::group_by(id_tcga_sample_barcode) %>%
    tidyr::spread(key=id_gene_symbol, value="quant_tmm_tpm") %>%
    tibble::column_to_rownames("id_tcga_sample_barcode") %>%
    cor() %>%
    .[lower.tri(.)] %>%
    mean()
  
  #fraction
  set_used_frac=length(intersect(unique(crc_tmm_gexp$id_gene_symbol),data_list[[i]]))/
    length(data_list[[i]])
  
  data.frame(id_geneset_name=i,quant_coherence_score=CS,set_used_frac=set_used_frac)
  
  
}
CS_scores_df

#### Mean signature score ####
# only take coherent signatures CS>0.2 & fraction of genes >0.5

CS_scores_df= CS_scores_df%>%
  dplyr::filter(quant_coherence_score > .2, set_used_frac > .5)

data_list=data_list[CS_scores_df$id_geneset_name]

mean_sig=foreach(i=names(data_list),.combine=rbind) %do% {
  
  crc_tmm_gexp %>%
    dplyr::filter(id_gene_symbol %in% data_list[[i]]) %>% 
    group_by(id_gene_symbol) %>%
    dplyr::mutate(quant_tmm_tpm = as.vector(scale(log2(quant_tmm_tpm + 1)))) %>%
    dplyr::group_by(id_tcga_sample_barcode) %>%
    dplyr::summarise(quant_signature_score=mean(quant_tmm_tpm)) %>%
    mutate(id_geneset_name=i)
}

# join coherence score and mean signature score
cs_mean_joined=left_join(mean_sig,CS_scores_df)

#### Signature scores from ROSETTASX ####
# done already by Julian Kreis 

# Fetch signature scores for signatures of RosettaSX with a high coherence score
# and enough available genes
gset = merck_tcga_rna_wts_gset() %>%
  dplyr::filter(quant_coherence_score > .2, set_used_frac > .5,
                id_tcga_cohort %in% c("BRCA"),
                meta_curator=="kreis") %>%
  dplyr::select(id_geneset_name,  id_tcga_sample_barcode, id_tcga_cohort,
                quant_signature_score, quant_coherence_score) %>%
  dplyr::collect() 

gset <- saveRDS(gset,"data/Tumor_bulk/gset.RDS")
gset <- readRDS("data/Tumor_bulk/gset.RDS")

gset$id_geneset_name=gsub("_\\d+$", "", gset$id_geneset_name)

# Remove internal gene sets e-g- Eike internal
gset=gset[!grepl("internal",gset$id_geneset_name), ]

# Join rosettasx signatures and newly added ones 
# Choose the column info to take 
col_take=c("id_tcga_sample_barcode",
           "id_geneset_name",
           "quant_signature_score")

all_df=rbind(cs_mean_joined[,col_take],
             gset[,col_take])

# Create the input matrix for heatmap
input_heat=
  acast(all_df, id_tcga_sample_barcode~id_geneset_name, value.var='quant_signature_score')


# Define color palettes
sig_palette = circlize::colorRamp2(breaks = c(-2, -.5, 0, .5, 2), 
                                   colors = c("#6473ed", "#2fd0e0", "#e8cbae",
                                              "#e6b83c", "#bf3255"))


#### Plotting ####

text_size=9
labelsize=10

# Plot also whole heatmap 
dim(input_heat)
cor_palette = circlize::colorRamp2(breaks = c(0.3, 0.5, 0.7), 
                                   colors = rev(c("#FFAA4C", "#ffffff", "#001E6C")))
IFNII_plot_all=rosettasx_heatmap_func(input_heat,
                                      covariance=T,
                                      correlation=F,
                                      selected_signature="IFN_II_Aybey",
                                      top=20, 
                                      textsize=text_size,
                                      label_size=labelsize,
                                      signature_colors=sig_palette,
                                      row_annot_colors=cor_palette)
IFNI_plot_all=rosettasx_heatmap_func(input_heat,
                                     covariance=T,
                                     correlation=F,
                                     selected_signature="IFN_I_Aybey",
                                     top=20, # top 20 correlated or covariance 
                                     textsize=text_size,
                                     # titlesize=text_size,
                                     label_size=labelsize,
                                     signature_colors=sig_palette,
                                     row_annot_colors=cor_palette)
# Print
pdf(paste0(save_dir,"TCGA_BRCA_top20_covariance_plots_1.pdf"),width = 6,height = 10)

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1)))  # 2 rows, 1 column

# Draw the first heatmap in the first row
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(IFNI_plot_all, newpage = FALSE)
popViewport()

# Draw the second heatmap in the second row
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
draw(IFNII_plot_all, newpage = FALSE)
popViewport(2)
dev.off()
