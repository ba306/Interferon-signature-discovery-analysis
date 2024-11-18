# Generate boxplots for gene expression profiles and gene expression signatures
# between responders and non-responders from 3 baseline ICB datasets
# Signatures: IFNI Aybey and IFNII Aybey

library(ProfilerAPI2)
library(dplyr)
library(foreach)
library(ggplot2)
library(rstatix)
library(ggpubr)
library(reshape2)
library(gridExtra)
library(grid)
library(readxl)

# Output directory
save_dir="Tumor_bulk_analysis/"
data_dir="data/Tumor_bulk/"

# Our IFN signatures
rosetta_new=readRDS(file = "published_signatures/IFN_signature_list_all_valid.RDS")

# Whole list comprimising genes and signatures for the downstream analysis
data_list=rosetta_new[names(rosetta_new) %in% "IFN_I_Aybey" |
                        names(rosetta_new) %in% "IFN_II_Aybey"    ]

# All genes to extract
all_genes=unique(unlist(data_list))

#### Datasets and analysis #### 

#### Imvigor210 Bladder ####
# Get expression data + metadata 

# get connection to profiler
api <- ProfilerAPI2::profiler_api(profile = "default")
# Get expression data from Profiler
imvigor <- 
  dplyr::tbl(api$conn,"genentech_imvigor210_rna_wts_gsnap_21q1_prod_jami_0") %>% 
  filter(gene_name %in% all_genes) %>%
  select(subject_id,gene_name,best_confirmed_overall_response,binaryresponse,gene_tpm,
         received_platinum,sample_collected_pre_platinum) %>%
  collect()

# remove samples without response info
imvigor=imvigor[!is.na(imvigor$binaryresponse),]

# Responders: SD/PD; non-responders: CR/PR 
imvigor$binaryresponse=
  ifelse(imvigor$binaryresponse=="SD/PD",
         "NR","R")

# Create boxplots for each gene or signature
# Cohort of platinum treated patients are chosen 
# Save individual plots as lists and compile together at the end

imvigor_plots=foreach(g=seq_along(data_list)) %do% {
  
  print(g)
  
  # Calculate mean signature scores (signatures) or gene expression (individual genes)
  # first log2(TPM+1) transformation and then average expression of every gene in a list
  
  score_Df=imvigor%>%
    filter(received_platinum=="Y") %>%
    filter(gene_name %in% unlist(data_list[g]))  %>%
    mutate(gene_logtpm = as.vector(scale(log2(1 + as.numeric(gene_tpm)))))%>%
    group_by(subject_id) %>%
    dplyr::summarize(score=mean(gene_logtpm, na.rm=TRUE),across())%>% 
    distinct(subject_id,.keep_all = T)
  
  # add p values between responders and non-responders
  stat.test <- score_Df %>%
    group_by(received_platinum) %>%
    t_test(score ~ binaryresponse) %>%
    # adjust_pvalue(method = "bonferroni") %>%
    add_significance()%>% add_xy_position(x = "binaryresponse")
  
  # Boxplot
  ggplot(score_Df, aes(x = binaryresponse, y = score)) +  
    geom_boxplot(outlier.shape = NA,color="black",width=0.3) + 
    geom_jitter(color="black", size=1, alpha=0.2) +theme_bw()+
    stat_pvalue_manual(stat.test)+ylab(names(data_list[g]))+
    xlab(NULL)+theme(text = element_text(size = 22))
  
}


# Other ICB datasets are
# downloaded from https://github.com/xmuyulab/ims_gene_signature

#### Kim gastric #####

# Expression data
kim_mtx=read.table(paste0(data_dir,"gas_korean_exp_data.csv"),
                   sep = ",",header = T,row.names = 1,check.names = F)
# Clinical data
kim_meta=read.table(paste0(data_dir,"gas_korean_cli_data.csv"),
                    sep = ",",header = T,row.names = 1)
# Binary response 
kim_meta$response=case_when(
  kim_meta$response=="-1" ~ "NR",
  kim_meta$response=="0" ~ "NR",
  kim_meta$response=="1" ~ "R",
  T~NA
)

# Select baseline samples
kim_meta=kim_meta%>%
  dplyr::filter(treatment=="pre")
# Select samples with exp and clinical data
kim_inter_pat=intersect(rownames(kim_meta),colnames(kim_mtx))
kim_meta=kim_meta[kim_inter_pat,]

# select baseline samples in expression data 
# z-scale the expression data
kim_mtx_scaled=scale(t(kim_mtx[,kim_inter_pat]))

# Calculate mean signature scores (signatures) or gene expression (individual genes)
# Add to the meta data

for(j in seq_along(data_list)){
  genes=intersect(data_list[[j]],colnames(kim_mtx_scaled))
  kim_meta[,names(data_list)[j]]= rowMeans(kim_mtx_scaled[,genes,drop=F],na.rm = T)
}

# Generate boxplots
kim_meta_plot_melt=melt(kim_meta)

kim_plots=foreach(g=seq_along(data_list)) %do% {
  print(g)
  
  # select gene or signature mean scores
  score_Df=kim_meta_plot_melt%>%
    filter(variable==names(data_list)[g])
  
  # add p values between responders and non-responders
  
  stat.test <- score_Df %>%
    t_test(value ~ response) %>%
    # adjust_pvalue(method = "bonferroni") %>%
    add_significance()%>% add_xy_position(x = "response")
  
  # Boxplot
  ggplot(score_Df, aes(x = response, y = value)) +  
    geom_boxplot(outlier.shape = NA,color="black",width=0.3) + 
    geom_jitter(color="black", size=1, alpha=0.2) +theme_bw()+
    stat_pvalue_manual(stat.test)+ylab(names(data_list[g]))+
    xlab(NULL)+theme(text = element_text(size = 22))
  
}

#### Van Allen melanoma ####
# Expression data
allen_mtx=read.table(paste0(data_dir,"mel_van_exp_data.csv"),
                     sep = ",",header = T,row.names = 1,check.names = F)
# Clinical data
allen_meta=read.table(paste0(data_dir,"mel_van_cli_data.csv"),
                      sep = ",",header = T,row.names = 1)
# Binary response 
allen_meta$response=case_when(
  allen_meta$response=="-1" ~ "NR",
  allen_meta$response=="0" ~ "NR",
  allen_meta$response=="1" ~ "R",
  T~NA
)

# Select baseline samples
allen_meta=allen_meta%>%
  dplyr::filter(treatment=="pre")
# Select samples with exp and clinical data
allen_inter_pat=intersect(rownames(allen_meta),colnames(allen_mtx))
allen_meta=allen_meta[allen_inter_pat,]

# select baseline samples in expression data 
# z-scale the expression data
allen_mtx_scaled=scale(t(allen_mtx[,allen_inter_pat]))

# Calculate mean signature scores (signatures) or gene expression (individual genes)
# Add to the meta data

for(j in seq_along(data_list)){
  genes=intersect(data_list[[j]],colnames(allen_mtx_scaled))
  allen_meta[,names(data_list)[j]]= rowMeans(allen_mtx_scaled[,genes,drop=F],na.rm = T)
}

# Generate boxplots
allen_meta_plot_melt=melt(allen_meta)

van_plots=foreach(g=seq_along(data_list)) %do% {
  print(g)
  
  # select gene or signature mean scores
  score_Df=allen_meta_plot_melt%>%
    filter(variable==names(data_list)[g])
  
  # add p values between responders and non-responders
  stat.test <- score_Df %>%
    t_test(value ~ response) %>%
    # adjust_pvalue(method = "bonferroni") %>%
    add_significance()%>% add_xy_position(x = "response")
  
  # Boxplot
  ggplot(score_Df, aes(x = response, y = value)) +  
    geom_boxplot(outlier.shape = NA,color="black",width=0.3) + 
    geom_jitter(color="black", size=1, alpha=0.2) +theme_bw()+
    stat_pvalue_manual(stat.test)+ylab(names(data_list[g]))+
    xlab(NULL)+theme(text = element_text(size = 22))
  
}

# Arrange plots
kim_arranged <- do.call("grid.arrange", c(kim_plots, ncol=1))
van_arranged <- do.call("grid.arrange", c(van_plots, ncol=1))
imvigor_arranged <- do.call("grid.arrange", c(imvigor_plots, ncol=1))


# Define the plots and their corresponding titles in a list
plot_list <- list(
  list(plot = imvigor_plots[[1]], title = "IMvigor210 (Bladder)"),
  list(plot = kim_plots[[1]], title = "Kim (Gastric)"),
  list(plot = van_plots[[1]], title = "Van Allen (Melanoma)")
  
)

# Create an empty list to store the arranged plots
arranged_plots <- list()

# Loop through the plot_list and arrange the plots with titles
for (item in plot_list) {
  arranged_plots[[length(arranged_plots) + 1]] <- grid.arrange(
    item$plot, ncol = 1,
    top = textGrob(item$title, gp = gpar(fontsize = 16, fontface = "bold"))
  )
}

# Combine the arranged plots into a single grid arrangement
final_arrangement <- do.call(grid.arrange, c(arranged_plots, ncol = 3))

# Print the final grid arrangement
ifnii_plots=list(kim_plots[[2]],
                 van_plots[[2]],
                 imvigor_plots[[2]])

ifnii_plots <- do.call(grid.arrange, c(ifnii_plots, ncol = 3))

pdf(paste0(save_dir,"ICB_IFNa_IFNg_Aybey_boxplots_new.pdf"),
    width = 8,height = 9)

grid.arrange(
  final_arrangement,ifnii_plots,
  bottom = textGrob("Response groups", 
                    gp = gpar(fontsize = 14, fontface = "bold"))
)
dev.off()



