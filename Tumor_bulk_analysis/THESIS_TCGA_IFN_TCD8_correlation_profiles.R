library(dplyr)
library(ggplot2)
library(foreach)
library(correlation)
library(gridExtra)
library(ggpubr)
library(xopggexpdatasets)
library(ProfilerAPI2)

save_dir="Tumor_bulk_analysis/"
data_dir= "data/Tumor_bulk/"

# load IFN signatures
rosetta_new=readRDS(file = "published_signatures/IFN_signature_list_all_valid.RDS")

# load T CD8 published signatures 
T_cd8_sign=readRDS(file = "published_signatures/T_CD8_published_signatures_list.RDS")

# Combine IFN and T CD8 signatures
rosetta_new=c(rosetta_new,T_cd8_sign)

###### Gene expression matrix from all 32 TCGA cohorts #####

# focus only on those genes we use for the analysis

all_genes=unique(unlist(rosetta_new))
length(all_genes) #532


# Load TCGA data from all cohorts only from genes to be analyzed
api <- ProfilerAPI2::profiler_api(profile = "default")

tcga <- getRefClass("TCGAExprDataset")$new(connection = api$conn)

data <- tcga$get_data_frame(c(all_genes,"meta_tcga_cohort","id_tcga_participant_barcode"))


# Calculate zscale mean signature scores for each cohort separately

mean_sig=foreach(j=unique(data$meta_tcga_cohort))%do%{
  print(j)
  cohort_means=foreach(i=seq_along(rosetta_new),.combine=cbind) %do% {
    intergenes=intersect(colnames(data),rosetta_new[[i]])
    a= log2(data[data$meta_tcga_cohort %in% j,intergenes]+1 ) %>% 
      apply( 2, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y)))%>% rowMeans()
    
  }
  colnames(cohort_means) =names(rosetta_new)
  cbind(cohort_means,data[data$meta_tcga_cohort %in% j,c("id_tcga_participant_barcode")])
}
names(mean_sig)=unique(data$meta_tcga_cohort)


# Our IFNs vs other IFN-Is 

IFN_I_names=names(rosetta_new)[!(grepl("IFNg",names(rosetta_new)) |
                                   grepl("T_CD8",names(rosetta_new)) |   
                                   grepl("Aybey",names(rosetta_new)) )
]

IFN_II_names=names(rosetta_new)[!(!grepl("IFNg",names(rosetta_new)) |
                                    grepl("T_CD8",names(rosetta_new)) |  
                                    grepl("IFNa",names(rosetta_new)) |
                                    grepl("IFNg_Ayers_extended",names(rosetta_new)) |  
                                    grepl("Aybey",names(rosetta_new)) )
]

T_CD8_names=names(rosetta_new)[grepl("T_CD8",names(rosetta_new))]

allIFNg=rosetta_new[!(!grepl("IFNg",names(rosetta_new)) |
                        grepl("T_CD8",names(rosetta_new)) |  
                        grepl("IFNa",names(rosetta_new)) |
                        grepl("Aybey",names(rosetta_new)) )
] %>% unlist() %>% as.vector() %>%unique()

intersect(allIFNg,rosetta_new$IFN_II_Aybey) #5 genes in common -.- GBP2 not common 

############ Correlation between IFNI or T CD8 and our IFN signature scores ######

# Calculate signature score correlations for each cohort
cor_df=foreach(j= unique(data$meta_tcga_cohort))%do% {
  print(j)
  mean_sig[[j]][,!grepl("id_tcga_participant_barcode",colnames(mean_sig[[j]]))] %>%na.omit()%>%
    cor(method = "pearson", use = "pairwise.complete.obs")
  
}
names(cor_df)=unique(data$meta_tcga_cohort)



#1. ifni vs ifnii: Can we separate IFNI and IFNII signal in TCGA using our signatures?
#2. tcd8 vs ifnii: Do we see similar distribution to TCD8 signatures?

df_ourIFN_vs_IFNI=foreach(i=names(cor_df),.combine=rbind) %do% {
  print(i)
  a=cor_df[[i]]
  
  #ifni--ifni
  b=a[IFN_I_names,IFN_I_names]
  b=b[lower.tri(b, diag = FALSE)]
  
  #ifnii--ifnii
  ifnii=a[IFN_II_names,IFN_II_names]
  ifnii=ifnii[lower.tri(ifnii, diag = FALSE)]
  
  d=rbind(b %>%as.vector() %>%data.frame(comparison="IFN-I- IFN-I signatures"),
          ifnii %>%as.vector() %>%data.frame(comparison="IFN-II- IFN-II signatures"),
          a[IFN_I_names,"IFN_I_Aybey"] %>%as.vector() %>%data.frame(comparison="IFN-I- IFN I Aybey"),
          a[IFN_I_names,"IFN_II_Aybey"] %>%as.vector() %>%data.frame(comparison="IFN-I- IFN II Aybey"),
          a[IFN_II_names,"IFN_I_Aybey"] %>%as.vector() %>%data.frame(comparison="IFN-II- IFN I Aybey"),
          a[IFN_II_names,"IFN_II_Aybey"] %>%as.vector() %>%data.frame(comparison="IFN-II- IFN II Aybey"),
          a["IFN_I_Aybey","IFN_II_Aybey"] %>%as.vector() %>%data.frame(comparison="IFN-I Aybey- IFN II Aybey"),
          a[IFN_II_names,IFN_I_names] %>%as.vector() %>%data.frame(comparison="IFN-I - IFN II")
          
          
          
          
  )
  
  
  # tcd8
  tcd8_df=a[T_CD8_names,T_CD8_names]
  tcd8_df=tcd8_df[lower.tri(tcd8_df, diag = FALSE)]
  
  tcd8_df_all=rbind(tcd8_df %>%as.vector() %>%data.frame(comparison="T CD8 - T CD8 signatures"),
                    a[T_CD8_names,"IFN_I_Aybey"] %>%as.vector() %>%data.frame(comparison="T CD8 - IFN I Aybey"),
                    a[T_CD8_names,"IFN_II_Aybey"] %>%as.vector() %>%data.frame(comparison="T CD8 - IFN II Aybey"),
                    a[T_CD8_names,IFN_II_names] %>%as.vector() %>%data.frame(comparison="T CD8 - IFN II")
                    
  )
  d_all=rbind(d,tcd8_df_all)
  data.frame(d_all,cohort=i)
  
  
}
colnames(df_ourIFN_vs_IFNI)[1]="cor"

# zCorrelation calculated 
df_ourIFN_vs_IFNI$z_cor=z_fisher(df_ourIFN_vs_IFNI$cor)

# t test function for comparisons

ttest_histo_p= function(data=filt_df){
  comp=unique(data$comparison)
  
  a=outer(
    comp, comp,
    Vectorize( function(g1,g2) {
      t.test(
        data$z_cor[ data$comparison == g1 ],
        data$z_cor[ data$comparison == g2 ]
      )$p.value
    } )
  )
  rownames(a)=comp
  colnames(a)=comp
  a}

unq_comp=c("IFN-I- IFN I Aybey" ,  "IFN-I- IFN II Aybey" )

# IFNI-IFNI vs IFNI-IFNI Aybey
# IFNI IFNI vs IFNI-IFNII Aybey
IFN_I_vs_ours_plot=foreach(i=unq_comp) %do% {
  filt_df=df_ourIFN_vs_IFNI %>% 
    filter(comparison %in%c("IFN-I- IFN-I signatures",i) )
  a=ttest_histo_p()
  
  d=mean(filt_df[filt_df$comparison==unique(filt_df$comparison)[1],"z_cor"])-
    mean(filt_df[filt_df$comparison==unique(filt_df$comparison)[2],"z_cor"])
  
  filt_df$comparison=as.factor(filt_df$comparison)
  filt_df$comparison=relevel(filt_df$comparison, ref = "IFN-I- IFN-I signatures")
  signig=  ifelse(a[1,2]<0.001,"***",
                  ifelse(a[1,2]<0.01,"**",
                         ifelse(a[1,2]<0.05,"*",""
                         )
                  )
  )
  
  gghistogram(filt_df, x = "z_cor", y = "..density..",
              add = "mean", rug = TRUE,
              fill = "comparison",
              add_density = TRUE)+
    annotate(geom="text", x=1.85, y=1.9, col="black", 
             label=paste0("d=",abs(round(d,2))))+
    annotate(geom="text", x=1.85, y=2, col="black",
             label=signig)
}


# IFNII-IFNII vs IFNII-IFNI Aybey
# IFNII IFNII vs IFNII-IFNII Aybey
unq_comp=c("IFN-II- IFN I Aybey" ,  "IFN-II- IFN II Aybey" )

IFN_II_vs_ours_plot=foreach(i=unq_comp) %do% {
  filt_df=df_ourIFN_vs_IFNI %>% 
    filter(comparison %in%c("IFN-II- IFN-II signatures",i) )
  a=ttest_histo_p()
  
  d=mean(filt_df[filt_df$comparison==unique(filt_df$comparison)[1],"z_cor"])-
    mean(filt_df[filt_df$comparison==unique(filt_df$comparison)[2],"z_cor"])
  
  filt_df$comparison=as.factor(filt_df$comparison)
  filt_df$comparison=relevel(filt_df$comparison, ref = "IFN-II- IFN-II signatures")
  signig=  ifelse(a[1,2]<0.001,"***",
                  ifelse(a[1,2]<0.01,"**",
                         ifelse(a[1,2]<0.05,"*",""
                         )
                  )
  )
  
  gghistogram(filt_df, x = "z_cor", y = "..density..",
              add = "mean", rug = TRUE,
              fill = "comparison",
              add_density = TRUE)+
    annotate(geom="text", x=1.85, y=1.9, col="black", 
             label=paste0("d=",abs(round(d,2))))+
    annotate(geom="text", x=1.85, y=2, col="black",
             label=signig)
}

IFN_II_vs_ours_plot[[1]]
IFN_II_vs_ours_plot[[2]]

# IFNI-IFNI aybey  vs IFNII-IFNI aybey

filt_df=df_ourIFN_vs_IFNI %>% 
  filter(comparison %in%c("IFN-I- IFN I Aybey","IFN-II- IFN I Aybey") )
a=ttest_histo_p()

d=mean(filt_df[filt_df$comparison==unique(filt_df$comparison)[1],"z_cor"])-
  mean(filt_df[filt_df$comparison==unique(filt_df$comparison)[2],"z_cor"])

signig=  ifelse(a[1,2]<0.001,"***",
                ifelse(a[1,2]<0.01,"**",
                       ifelse(a[1,2]<0.05,"*",""
                       )
                )
)

gghistogram(filt_df, x = "z_cor", y = "..density..",
            add = "mean", rug = TRUE,
            fill = "comparison",
            add_density = TRUE)+
  annotate(geom="text", x=1.85, y=1.9, col="black", 
           label=paste0("d=",abs(round(d,2))))+
  annotate(geom="text", x=1.85, y=2, col="black",
           label=signig)

# IFNI-IFNII aybey vs IFNII-IFNII aybey

filt_df=df_ourIFN_vs_IFNI %>% 
  filter(comparison %in%c("IFN-I- IFN II Aybey","IFN-II- IFN II Aybey") )
a=ttest_histo_p()

d=mean(filt_df[filt_df$comparison==unique(filt_df$comparison)[1],"z_cor"])-
  mean(filt_df[filt_df$comparison==unique(filt_df$comparison)[2],"z_cor"])

signig=  ifelse(a[1,2]<0.001,"***",
                ifelse(a[1,2]<0.01,"**",
                       ifelse(a[1,2]<0.05,"*",""
                       )
                )
)

gghistogram(filt_df, x = "z_cor", y = "..density..",
            add = "mean", rug = TRUE,
            fill = "comparison",
            add_density = TRUE)+
  annotate(geom="text", x=1.85, y=1.9, col="black", 
           label=paste0("d=",abs(round(d,2))))+
  annotate(geom="text", x=1.85, y=2, col="black",
           label=signig)

# IFNI aybey IFNII aybey vs IFNI-IFNII
filt_df=df_ourIFN_vs_IFNI %>% 
  filter(comparison %in%c("IFN-I Aybey- IFN II Aybey","IFN-I - IFN II") )
a=ttest_histo_p()

d=mean(filt_df[filt_df$comparison==unique(filt_df$comparison)[1],"z_cor"])-
  mean(filt_df[filt_df$comparison==unique(filt_df$comparison)[2],"z_cor"])

signig=  ifelse(a[1,2]<0.001,"***",
                ifelse(a[1,2]<0.01,"**",
                       ifelse(a[1,2]<0.05,"*",""
                       )
                )
)

gghistogram(filt_df, x = "z_cor", y = "..density..",
            add = "mean", rug = TRUE,
            fill = "comparison",
            add_density = TRUE)+
  annotate(geom="text", x=1.85, y=1.9, col="black", 
           label=paste0("d=",abs(round(d,2))))+
  annotate(geom="text", x=1.85, y=2, col="black",
           label=signig)


# T CD8 analysis

unq_comp=unique(df_ourIFN_vs_IFNI$comparison)
unq_comp=unq_comp[grepl("T CD8",unq_comp)][-1]

cd8_vs_ours_plot=foreach(i=unq_comp) %do% {
  filt_df=df_ourIFN_vs_IFNI %>% 
    filter(comparison %in%c("T CD8 - T CD8 signatures",i) )
  a=ttest_histo_p()
  
  d=mean(filt_df[filt_df$comparison==unique(filt_df$comparison)[1],"z_cor"])-
    mean(filt_df[filt_df$comparison==unique(filt_df$comparison)[2],"z_cor"])
  
  filt_df$comparison=as.factor(filt_df$comparison)
  filt_df$comparison=relevel(filt_df$comparison, ref = "T CD8 - T CD8 signatures")
  signig=  ifelse(a[1,2]<0.001,"***",
                  ifelse(a[1,2]<0.01,"**",
                         ifelse(a[1,2]<0.05,"*","ns"
                         )
                  )
  )
  
  gghistogram(filt_df, x = "z_cor", y = "..density..",
              add = "mean", rug = TRUE,
              fill = "comparison",
              add_density = TRUE)+
    annotate(geom="text", x=1.85, y=1.9, col="black", 
             label=paste0("d=",round(d,2)))+
    annotate(geom="text", x=1.85, y=2, col="black",
             label=signig)
}



# Arrange the first set of plots in a grid
plot1 <- grid.arrange(IFN_I_vs_ours_plot[[1]],
                      IFN_II_vs_ours_plot[[2]],
                      IFN_I_vs_ours_plot[[2]],
                      IFN_II_vs_ours_plot[[1]], 
                      ncol = 2)

# Arrange the second set of plots in a grid
plot2 <- grid.arrange(cd8_vs_ours_plot[[1]], cd8_vs_ours_plot[[2]], ncol = 2)

# Combine the two grids into one
final_plot <- ggarrange(plot1, plot2, 
                        heights = c(2, 1), 
                        labels = c("A", "B"), 
                        font.label = list(size = 24),ncol = 1)

pdf(paste0(save_dir,"zCorhistos_IFNI_cd8_vs_ourIFNs.pdf"),width = 25,height = 20)

final_plot

dev.off()

# saved from plots window zCovariance.pdf better resolution"

