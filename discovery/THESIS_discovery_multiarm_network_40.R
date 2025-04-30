#IFN I vs IFN II signature discovery workflow
# Based on network meta analysis approach
# IFN IFN comparisons added 
# 40 -makecluster

library(limma)
library(netmeta)
library(DESeq2)
library(dplyr)
library(foreach)
library(doParallel)
library(stringr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(gridtext)
library(grid)
library(gplots)

save_dir="./discovery/"
data_dir="./data/IFN_bulk/" 

# Load discovery datasets incl. metadata and expression data
load(paste0(data_dir,
            "discovery_exp_meta_rdata.RData"))

# Operate only on common genes
inter_5=Reduce(intersect, list(rownames(ziegler_c),
                               rownames(jank_c),
                               rownames(colli_c),
                               rownames(rai_c),
                               rownames(fuji_c)
))
length(inter_5)
#11484

# Choose only those common genes
ziegler_g=ziegler_c[inter_5,]
jank_g=jank_c[inter_5,]
fuji_g=fuji_c[inter_5,]
colli_g=colli_c[inter_5,]
rai_g=rai_c[inter_5,]


#### Rank genes by their variations and choose top 5k common most HVGs ####
all_list=list(
  ziegler_g,
  jank_g,
  colli_g,
  rai_g,
  fuji_g)

# rank by variance in each dataset 
all_list_var=lapply(all_list, function(x){
  x=as.data.frame(x)
  x$var=apply(x, 1, var)
  x$rank_var=rank(x$var)
  x[,"rank_var",drop=F]
})

# sum of variance rankings
var_rank_sum=Reduce(all_list_var,f = "+")

# choose 5k highest
topvar5k=var_rank_sum%>%
  dplyr::arrange_(~ desc(rank_var)) %>%
  dplyr::slice(1:5000)

ziegler_g=ziegler_g[rownames(topvar5k),]
jank_g=jank_g[rownames(topvar5k),]
fuji_g=fuji_g[rownames(topvar5k),]
colli_g=colli_g[rownames(topvar5k),]
rai_g=rai_g[rownames(topvar5k),]

#### Differential gene analysis ####
### study1 Jankowski : treatments versus control
group = jank_df$Stim
group=as.factor(group)
group=relevel(group,"untreated")

M = model.matrix(~ group)
fit = lmFit(jank_g, M)
fit = eBayes(fit,trend = T)

jank_list=foreach(i=2:ncol(fit$p.value)) %do% {
  p.S1_a = fit$p.value[,i]
  fc.S1_a = fit$coefficients[,i]
  fce.S1_a = sqrt(fit$s2.post) * sqrt(fit$cov.coefficients[i,i])
  data.frame(fc=fc.S1_a,fc_se=fce.S1_a,pval=p.S1_a)
}
names(jank_list)=levels(group)[-1]

group = jank_df$Stim
group=as.factor(group)
group=relevel(group,"IFNa")

M = model.matrix(~ group)
fit = lmFit(jank_g, M)
fit = eBayes(fit,trend = T)

jank_list_IFNa=foreach(i=2:ncol(fit$p.value)) %do% {
  p.S1_a = fit$p.value[,i]
  fc.S1_a = fit$coefficients[,i]
  fce.S1_a = sqrt(fit$s2.post) * sqrt(fit$cov.coefficients[i,i])
  data.frame(fc=fc.S1_a,fc_se=fce.S1_a,pval=p.S1_a)
}
names(jank_list_IFNa)=levels(group)[-1]

fit$stdev.unscaled %>% head()
jank_list_IFNa=jank_list_IFNa[names(jank_list_IFNa)[names(jank_list_IFNa)!="untreated"]]
names(jank_list_IFNa)


group = jank_df$Stim
group=as.factor(group)
group=relevel(group,"IFNb")

M = model.matrix(~ group)
fit = lmFit(jank_g, M)
fit = eBayes(fit,trend = T)

jank_list_IFNb=foreach(i=2:ncol(fit$p.value)) %do% {
  p.S1_a = fit$p.value[,i]
  fc.S1_a = fit$coefficients[,i]
  fce.S1_a = sqrt(fit$s2.post) * sqrt(fit$cov.coefficients[i,i])
  data.frame(fc=fc.S1_a,fc_se=fce.S1_a,pval=p.S1_a)
}
names(jank_list_IFNb)=levels(group)[-1]

jank_list_IFNb=jank_list_IFNb[names(jank_list_IFNb)[names(jank_list_IFNb)!="untreated"]]
jank_list_IFNb=jank_list_IFNb[names(jank_list_IFNb)[names(jank_list_IFNb)!="IFNa"]]

names(jank_list_IFNb)



### study2 Ziegler Donor 1: treatment B versus control
group = ziegler_df$Stim
group=as.factor(group)
group=relevel(group,"untreated")


M = model.matrix(~ group)
fit = lmFit(ziegler_g, M)
fit = eBayes(fit,trend = T)

ziegler_list=foreach(i=2:ncol(fit$p.value)) %do% {
  print(i)
  p.S1_a = fit$p.value[,i]
  fc.S1_a = fit$coefficients[,i]
  fce.S1_a = sqrt(fit$s2.post) * sqrt(fit$cov.coefficients[i,i])
  data.frame(fc=fc.S1_a,fc_se=fce.S1_a,pval=p.S1_a)
}
names(ziegler_list)=levels(group)[-1]

lapply(ziegler_list, function(x){summary(x$fc)})


group=ziegler_df %>%
  filter(Stim!="untreated") 
pat=rownames(group)
group=as.factor(group$Stim)
group=relevel(group,"IFNa")

Dose=ziegler_df %>%
  filter(Stim!="untreated") 
Dose =Dose$Dose
Dose=as.factor(Dose)
# group=relevel(group,"untreated")


M = model.matrix(~ group+Dose)

fit = lmFit(ziegler_g[,pat], M)
fit = eBayes(fit,trend = T)

ziegler_list_IFNa=foreach(i=2) %do% {
  p.S1_a = fit$p.value[,i]
  fc.S1_a = fit$coefficients[,i]
  fce.S1_a = sqrt(fit$s2.post) * sqrt(fit$cov.coefficients[i,i])
  data.frame(fc=fc.S1_a,fc_se=fce.S1_a,pval=p.S1_a)
}
names(ziegler_list_IFNa)=levels(group)[-1]


#study 3 Fujiwara
group = fuji_df$Stim
group=as.factor(group)
group=relevel(group,"untreated")

M = model.matrix(~ group)
fit = lmFit(fuji_g, M)
fit = eBayes(fit,trend = T)

fuji_list=foreach(i=2:ncol(fit$p.value)) %do% {
  print(i)
  p.S1_a = fit$p.value[,i]
  fc.S1_a = fit$coefficients[,i]
  fce.S1_a = sqrt(fit$s2.post) * sqrt(fit$cov.coefficients[i,i])
  dt=data.frame(fc=fc.S1_a,fc_se=fce.S1_a,pval=p.S1_a)
  rownames(dt)=names(fce.S1_a)
  dt
}
names(fuji_list)=levels(group)[-1]


#study 4 Colli
group = colli_df$Stim
group=as.factor(group)
group=relevel(group,"untreated")

time=colli_df$time
time=as.factor(time)
time <- factor(time, levels = c(2,8,18))

M = model.matrix(~ group+time)
fit = lmFit(colli_g, M)
fit = eBayes(fit,trend = T)

colli_list=foreach(i=2) %do% {
  print(i)
  p.S1_a = fit$p.value[,i]
  fc.S1_a = fit$coefficients[,i]
  fce.S1_a = sqrt(fit$s2.post) * sqrt(fit$cov.coefficients[i,i])
  dt=data.frame(fc=fc.S1_a,fc_se=fce.S1_a,pval=p.S1_a)
  rownames(dt)=names(p.S1_a)
  dt
}
names(colli_list)=levels(group)[-1]

#study 5 Rai
group = rai_df$Stim
group=as.factor(group)
group=relevel(group,"untreated")

M = model.matrix(~ group)
fit = lmFit(rai_g, M)
fit = eBayes(fit,trend = T)

rai_list=foreach(i=2:ncol(fit$p.value)) %do% {
  print(i)
  p.S1_a = fit$p.value[,i]
  fc.S1_a = fit$coefficients[,i]
  fce.S1_a = sqrt(fit$s2.post) * sqrt(fit$cov.coefficients[i,i])
  dt=data.frame(fc=fc.S1_a,fc_se=fce.S1_a,pval=p.S1_a)
  rownames(dt)=names(fce.S1_a)
  dt
}
names(rai_list)=levels(group)[-1]


fit_list = c(jank_list,
             jank_list_IFNa,
             jank_list_IFNb,
             ziegler_list,
             ziegler_list_IFNa,
             fuji_list,
             colli_list,
             rai_list)

names(fit_list)=c(paste0("jankowski_",names(jank_list)),
                  paste0("jankowski_IFNa_vs_",names(jank_list_IFNa)),
                  paste0("jankowski_IFNb_vs_",names(jank_list_IFNb)),
                  
                  paste0("ziegler_",names(ziegler_list)),
                  paste0("ziegler_IFNa_vs_",names(ziegler_list_IFNa)),
                  paste0("fujiwara_",names(fuji_list)),
                  paste0("colli_",names(colli_list)),
                  paste0("rai_",names(rai_list)))


# P value histograms between different conditions in each discovery dataset

hist_list=foreach(i=seq_along(fit_list)) %do% {
  rank_test=fit_list[[i]]
  
  gname=if(!grepl("vs",names(fit_list)[i])){
    paste0(names(fit_list)[i],"_vs_untreated")
  }else{
    names(fit_list)[i]
  }
  gname= gsub("_"," ",gname)

  gname=str_to_title(gname) 
  gname= gsub("Vs","vs",gname)
  
  gname= gsub("Ifn","IFN",gname)
  
  
  ggplot(rank_test, aes(x=pval)) +
    geom_histogram(bins = 20,fill="white", color="black")+
    #scale_color_grey()+#scale_fill_grey() +
    theme_classic()+
    # geom_hline(yintercept =nrow(rank_test)/20 ,color="red",
    #            linetype="dashed")+
    ggtitle(
                   gname)+ 
    theme(plot.title = element_text(size = 20),
          axis.text=element_text(size=18),
          axis.title=element_text(size=18))
}

# Axis labels
yleft = richtext_grob("Count", rot=90,gp= gpar(fontsize = 22))

bottom = richtext_grob(
  text = "p-value",gp= gpar(fontsize = 22)
)

# Lay out plots
hist_final_plot=grid.arrange(grobs=hist_list%>% map(~.x + labs(x=NULL, y=NULL)), ncol =3, 
                    left = yleft, bottom = bottom)

ggsave(paste0(save_dir2,"p_val_distribution_disc_datasets.pdf"), width = 14,height = 12,scale=1,
       plot = hist_final_plot)

#### Network meta-analysis ####

# Treatment groups 
treat2_all = c(paste0("Untreated_vs_",names(jank_list)),
               paste0("IFNa_vs_",names(jank_list_IFNa)),
               paste0("IFNb_vs_",names(jank_list_IFNb)),
               
               paste0("Untreated_vs_",names(ziegler_list)),
               paste0("IFNa_vs_",names(ziegler_list_IFNa)),
               paste0("Untreated_vs_",names(fuji_list)),
               paste0("Untreated_vs_",names(colli_list)),
               paste0("Untreated_vs_",names(rai_list)))


treat2=gsub("(.*)_vs_(.*)","\\2",treat2_all)

treat1=gsub("(.*)_vs_(.*)","\\1",treat2_all)

# study groups
studlab = c(rep("jank",length(jank_list)),
            rep("jank",length(jank_list_IFNa)),
            rep("jank",length(jank_list_IFNb)),
            
            rep("ziegler",length(ziegler_list)),
            rep("ziegler",length(ziegler_list_IFNa)),
            rep("fuji",length(fuji_list)),
            rep("colli",length(colli_list)),
            rep("rai",length(rai_list)))

# Faster implementation using parallel computing

cl <- parallel::makeCluster(40)
doParallel::registerDoParallel(cl)
num=rownames(jank_list$IFNb)

# IFN gene lists comparing different treatments
IFN_list=list('a_b'=list(),"a_g"=list(),"b_g"=list()) #for treatment vs treatment 
IFN_control_list=list('a'=list(),"b"=list(),"g"=list())  # for treatment vs control 

print("started")

All_network_list=foreach (j=seq_along(num)) %do%  {
  print(j)
  d1=dplyr::bind_rows(lapply(jank_list,function(x)x[num[j],-3]))
  d2=dplyr::bind_rows(lapply(jank_list_IFNa,function(x)x[num[j],-3]))
  d3=dplyr::bind_rows(lapply(jank_list_IFNb,function(x)x[num[j],-3]))
  
  d4=dplyr::bind_rows(lapply(ziegler_list,function(x)x[num[j],-3]))
  d5=dplyr::bind_rows(lapply(ziegler_list_IFNa,function(x)x[num[j],-3]))
  
  d6=dplyr::bind_rows(lapply(fuji_list,function(x)x[num[j],-3]))
  d7=dplyr::bind_rows(lapply(colli_list,function(x)x[num[j],-3]))
  d8=dplyr::bind_rows(lapply(rai_list,function(x)x[num[j],-3]))
  
  ds=rbind(d1,d2,d3,d4,d5,d6,d7,d8)
  
  ds$treat1=treat1
  ds$treat2=treat2
  ds$studlab=studlab
  gene_name=num[j]
  
  net1 = netmeta::netmeta(fc, fc_se, treat1, treat2, studlab, data=ds, sm="MD",reference.group = "Untreated")
  S =  summary(net1)

  # Treatment vs treatment  lists 
  
  IFN_list[[1]]=rbind(IFN_list[[1]],
                      data.frame(logfc=S$random$TE["IFNb","IFNa"],
                                 pval=S$random$p["IFNb","IFNa"],
                                 row.names = gene_name))
  IFN_list[[2]]=rbind(IFN_list[[2]],
                      data.frame(logfc=S$random$TE["IFNg","IFNa"],
                                 pval=S$random$p["IFNg","IFNa"],
                                 row.names = gene_name))
  IFN_list[[3]]=rbind(IFN_list[[3]],
                      data.frame(logfc=S$random$TE["IFNg","IFNb"],
                                 pval=S$random$p["IFNg","IFNb"],
                                 row.names = gene_name))
  # Treatment vs  control lists 
  IFN_control_list[[1]]=rbind(IFN_control_list[[1]],
                              data.frame(logfc=S$random$TE["Untreated","IFNa"],
                                         pval=S$random$p["Untreated","IFNa"],
                                         row.names = gene_name))
  IFN_control_list[[2]]=rbind(IFN_control_list[[2]],
                              data.frame(logfc=S$random$TE["Untreated","IFNb"],
                                         pval=S$random$p["Untreated","IFNb"],
                                         row.names = gene_name))
  IFN_control_list[[3]]=rbind(IFN_control_list[[3]],
                              data.frame(logfc=S$random$TE["Untreated","IFNg"],
                                         pval=S$random$p["Untreated","IFNg"],
                                         row.names = gene_name))
  
  net1
  
}
parallel::stopCluster(cl)

print("finished")

# save two list separately 
save(IFN_list,file = paste0(save_dir,"IFN_list_treatment_interaction.RDS"))
save(IFN_control_list,file = paste0(save_dir,"IFN_list_control_interaction.RDS"))

load(file = paste0(save_dir,"IFN_list_treatment_interaction.RDS"))
load(file = paste0(save_dir,"IFN_list_control_interaction.RDS"))


#### Filter genes ####
# Based on p values and fold change filter genes 

IFN_list_rank=lapply(IFN_list,function(x){
  
  x=x[abs(x$logfc)>2.5,]
  x$padj=p.adjust(x$pval,method = "BH")
  x=x[x$padj<0.05,]
  
  x$gene=rownames(x)
  x
})


v.table <- venn(lapply(IFN_list_rank, function(x){x$gene}))
lapply(IFN_list_rank, nrow)
# $a_b
# [1] 16
# 
# $a_g
# [1] 26
# 
# $b_g
# [1] 79

IFN_list_contol_rank=lapply(IFN_control_list,function(x){
  
  x=x[x$logfc>3,]
  x$padj=p.adjust(x$pval,method = "BH")
  x=x[x$padj<0.05,]
  
  x$gene=rownames(x)
  x
})


lapply(IFN_list_contol_rank, nrow)
v.table <- venn(lapply(IFN_list_contol_rank, function(x){x$gene}))
# $a
# [1] 55
# 
# $b
# [1] 111
# 
# $g
# [1] 38

# saveRDS(IFN_list_rank,file = paste0(save_dir,"IFN_list_filtered_interaction.RDS"))

# Final signatures 
# Positive log fold change genes for each comparison

# IFN_list_rank=readRDS(paste0(save_dir,"IFN_list_filtered_interaction.RDS"))
IFN_list_ind=list()
IFN_list_ind$a=c(IFN_list_rank$a_b[IFN_list_rank$a_b$logfc>0,"gene"],
                 IFN_list_rank$a_g[IFN_list_rank$a_g$logfc>0,"gene"])
IFN_list_ind$b=c(IFN_list_rank$a_b[IFN_list_rank$a_b$logfc<0,"gene"],
                 IFN_list_rank$b_g[IFN_list_rank$b_g$logfc>0,"gene"])
IFN_list_ind$g=c(IFN_list_rank$a_g[IFN_list_rank$a_g$logfc<0,"gene"],
                 IFN_list_rank$b_g[IFN_list_rank$b_g$logfc<0,"gene"])

lapply(IFN_list_ind, length)
# $a
# [1] 20
# 
# $b
# [1] 94
# 
# $g
# [1] 7

v.table <- venn(IFN_list_ind)

# intersected genes between different IFN lists
for(i in 1:length(IFN_list_rank)){
  # print(IFN_list_ind[[i]])
  IFN_list_ind[[i]]=intersect(IFN_list_ind[[i]],rownames(IFN_list_contol_rank[[i]])) %>% print()
  
}

lapply(IFN_list_ind, length)
v.table <- venn(IFN_list_ind)

# $a
# [1] 20
# 
# $b
# [1] 60
# 
# $g
# [1] 6


# Unique genes for each comparison 

IFN_list_ind$a=unique(IFN_list_ind$a)
IFN_list_ind$b=unique(IFN_list_ind$b)
IFN_list_ind$g=unique(IFN_list_ind$g)

lapply(IFN_list_ind, length)
# [1] 20
# 
# $b
# [1] 60
# 
# $g
# [1] 6

# Intersection between different groups?
intersect(IFN_list_ind$a,IFN_list_ind$b) %>% length() #20 ----> all IFNa signatures in IFNb 
intersect(IFN_list_ind$b,IFN_list_ind$g) %>% length() #0 

# all IFNa genes are also present in IFNb so for IFN b remove those common genes
# IFNa signatures is here an IFN-I signature !!! 

a_b_common=intersect(IFN_list_ind$b,IFN_list_ind$a)
IFN_list_ind$b=IFN_list_ind$b[-which(IFN_list_ind$b %in%a_b_common )]

names(IFN_list_ind)=c("IFNa_Aybey","IFNb_Aybey","IFNg_Aybey")
IFN_list_ind_df=melt(IFN_list_ind)

write.table(IFN_list_ind_df,paste0(save_dir,"IFN_Aybey_signatures.txt"),sep="\t")
saveRDS(IFN_list_ind,file = paste(save_dir,"IFN_Aybey_finallist.RDS"))


lapply(IFN_list_ind, length)
# $IFNa_Aybey
# [1] 20
# 
# $IFNb_Aybey
# [1] 40
# 
# $IFNg_Aybey
# [1] 6

# Generate plots for boxplots and heatmap in discovery and validation datasets using IFN Aybey signatures and published ones
source(paste0(save_dir,"THESIS_boxplots_discovery_validation_wstats.R"))
source(paste0(save_dir,"THESIS_Discovery_heatmap_woIFNb.R"))
source(paste0(save_dir,"THESIS_Validation_heatmap_woIFNb.R"))

