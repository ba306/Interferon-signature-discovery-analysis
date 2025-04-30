library(foreach)
library(dplyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(gridtext)
source( "functions/CScores.R")

# Heatmaps for our IFN signatures along with published IFN signatures
# remove IFNb
# In discovery datasets

save_dir="discovery/" 
data_dir="data/IFN_bulk/" 

### Load IFN gene signatures

rosetta_new=readRDS(file ="published_signatures/IFN_signature_list_all_valid.RDS")

## Colors for stimulations 

colors=sapply(c("#e41a1c", #red
                "#377eb8",# blue
                "#4daf4a", #green
                "#984ea3", #purple
                "#ff7f00", #orange
                "#ffff33", #yellow
                "#a65628"), toupper) # brown
colors=c(colors[1],
         colors[2],
         colors[3],
         colors[4],
         colors[5])
names(colors)=c("IFNa","IFNb","IFNg","Control","IFNlambda")

text_size=50
axis_size=30
grid_width=1.5
annot_height=1.5

# Load discovery data  

load(
  file = paste0(data_dir,
                "discovery_exp_meta_rdata.RData"))

### Coherence score ####


#for fast calculation only genes from our gene lists used

gene_coherence=unlist(rosetta_new,use.names = F)

exp_list=list(
  ziegler_c,
  jank_c,
  colli_c,
  rai_c,
  fuji_c)
cs_list_ind= foreach(j=seq_along(exp_list),.combine = cbind) %do% {
  print(j)
  foreach(i=seq_along(rosetta_new),.combine = rbind) %do% {
    print(i)
    gen= intersect(rownames(exp_list[[j]]),rosetta_new[[i]])
    pear=cor(t( exp_list[[j]][gen,]), method = "pearson")
    score=coherence_score_val(gen, pear)
    round(score,digits = 2)
  }
} 
rownames(cs_list_ind)=names(rosetta_new)
colnames(cs_list_ind)= c("ziegler",
                         "jankowski",
                         "colli",
                         "rai",
                         "fujiwara"
)

### Individual heatmaps ####

# Scale and calculate mean signature and also save corresponding
# meta (e.g. ziegler_df) file as meta 
mean_sig_function= function(data) {
  data=as.matrix(data)
  data= t(scale(t(data), center=TRUE, scale=TRUE))
  
  foreach(i=seq_along(names(rosetta_new)),.combine = cbind) %do%{
    print(i)
    inter_genes=intersect(rosetta_new[[i]],rownames(data))
    data[inter_genes,,drop=F] %>%
      colMeans(na.rm = T) %>%
      as.data.frame() %>%
      setNames(nm = names(rosetta_new)[i])
  }
}

# untreated changed to Control in metadata 

var_names <- names(as.list(.GlobalEnv))
df_names=var_names[grepl("(.)_df",var_names)]
for(i in df_names){
  var= get(i)
  var$Stim=gsub("untreated","Control",var$Stim)
  assign(x =i,var)
}

### Ziegler donor 1 #####

mean_sig= mean_sig_function(ziegler_c)

meta=ziegler_df
meta$Dose=as.numeric(meta$Dose)
meta=meta[order(meta$Stim),]
meta=meta[order(meta$Dose),]
meta=meta[order(meta$Stim),]

colors_select=colors[unique(meta$Stim)]

dose_col=rev(brewer.pal(n=7, name="RdBu"))
names(dose_col)=sort(unique(meta$Dose))


col_annot=HeatmapAnnotation(df = meta[,c("Stim","Dose"),drop=F],
                            annotation_legend_param = list(
                              Stim = list(direction = "horizontal",
                                          nrow=1, 
                                          labels_gp = gpar(fontsize = text_size),
                                          title_gp = gpar(fontsize = text_size, fontface = "bold"),
                                          grid_width = unit(grid_width, "cm")),
                              Dose=list(direction = "horizontal",
                                        nrow=1, 
                                        labels_gp = gpar(fontsize = text_size),
                                        title_gp = gpar(fontsize = text_size, fontface = "bold"),
                                        grid_width = unit(grid_width, "cm")
                                        )
                              ),
                            col=list(Stim=colors,
                                      Dose=dose_col),
                            annotation_name_gp = gpar(fontsize = text_size,fontface = "bold"),
                            simple_anno_size= unit(annot_height, "cm"))
col_annot@anno_list[["Dose"]]@color_mapping@levels=names(dose_col)
col_annot@anno_list[["Dose"]]@color_mapping@colors=dose_col

mean_scale= mean_sig[rownames(meta),]

row_ha = rowAnnotation(Coherence_score = anno_barplot(cs_list_ind[,"ziegler"], border = F,ylim = c(0,1),
                                                      axis_param = list(gp = gpar(fontsize = axis_size))  # Bold and larger font for axis labels
),
annotation_name_gp = gpar(fontsize = text_size, fontface = "bold"),
width = unit(2, "cm"))  # Bold and larger font for annotation name



ziegler_beas_p=Heatmap(t(mean_scale),cluster_columns = T , right_annotation = row_ha,
                       cluster_rows = T,name = "Mean signature score",
                       top_annotation= col_annot,column_title = "Ziegler donor1",
                       column_title_gp  = gpar(fontsize = text_size, fontface = "bold"),
                       show_column_names = F,
                       row_names_side = "right",row_dend_side = "left",
                       heatmap_legend_param = list(direction = "horizontal",
                                                   legend_width = unit(20, "cm"),
                                                   legend_gp  =gpar(fontsize = text_size),
                                                   labels_gp  = gpar(fontsize = text_size),
                                                   title_gp = gpar( fontsize = text_size)),
                       row_names_max_width = max_text_width(
                         colnames(mean_scale), 
                         gp = gpar(fontsize=text_size)
                       ),
                       row_names_gp = gpar(col = c(rep("red",2), 
                                                   rep("black",ncol(mean_scale)-2)),
                                           fontsize=text_size),
                       row_dend_width = unit(3, "cm"),
                       height = unit(40, "cm"),
                       width = unit(20, "cm")
)


draw(ziegler_beas_p, heatmap_legend_side = "top",
     annotation_legend_side = "bottom",newpage = T)

#### Jankowski ####
mean_sig= mean_sig_function(jank_c)

meta=jank_df
meta$Stim=meta$Stim
meta=meta[order(meta$Stim),]


colors_select=colors[unique(meta$Stim)]

col_annot=HeatmapAnnotation(df = meta[,c("Stim"),drop=F],
                            annotation_legend_param = list(
                              Stim = list(direction = "horizontal",
                                          grid_width = unit(grid_width, "cm"),
                              nrow=1, labels_gp = gpar(fontsize = text_size),
                              title_gp = gpar(fontsize = text_size, fontface = "bold")))
                            ,col=list(Stim=colors),
                            annotation_name_gp = gpar(fontsize = text_size,fontface = "bold"),
                            simple_anno_size= unit(annot_height, "cm"))

mean_scale= mean_sig[rownames(meta),]

row_ha = rowAnnotation(Coherence_score = anno_barplot(cs_list_ind[,"jankowski"], border = F,ylim = c(0,1),
                                                      axis_param = list(gp = gpar(fontsize = axis_size))  # Bold and larger font for axis labels
),
annotation_name_gp = gpar(fontsize = text_size, fontface = "bold"),
width = unit(2, "cm"))  # Bold and larger font for annotation name

jank_p=Heatmap(t(mean_scale),cluster_columns = T ,right_annotation = row_ha,
               cluster_rows = T,name = "Mean signature score",
               top_annotation= col_annot,column_title = "Jankowski",
               column_title_gp  = gpar(fontsize = text_size, fontface = "bold"),
               show_column_names = F,
               row_names_side = "right",row_dend_side = "left",
               heatmap_legend_param = list(direction = "horizontal",
                                           legend_width = unit(20, "cm"),
                                           legend_gp  =gpar(fontsize = text_size),
                                           labels_gp  = gpar(fontsize = text_size),
                                           title_gp = gpar( fontsize = text_size)),
               row_names_max_width = max_text_width(
                 colnames(mean_scale), 
                 gp = gpar(fontsize=text_size)
               ),
               row_names_gp = gpar(col = c(rep("red",2), 
                                           rep("black",ncol(mean_scale)-2)),
                                   fontsize=text_size),
               width = unit(10, "cm"),
               row_dend_width = unit(3, "cm"),
               height = unit(40, "cm")
)

draw(jank_p, heatmap_legend_side = "top",
     annotation_legend_side = "bottom",newpage = T)


##### Colli ####
mean_sig= mean_sig_function(colli_c)

meta=colli_df
meta$time=meta$time
meta$time=as.numeric(meta$time)

meta=meta[order(meta$Stim),]


colors_select=colors[unique(meta$Stim)]


dose_col=rev(brewer.pal(n=3, name="RdBu"))
names(dose_col)=sort(unique(meta$time))

colnames(meta)=gsub("time","Time",colnames(meta))

col_annot=HeatmapAnnotation(df = meta[,c("Stim","Time"),drop=F],
                            annotation_legend_param = list(
                              Stim = list(direction = "horizontal",
                                          nrow=1, 
                                          labels_gp = gpar(fontsize = text_size),
                                          title_gp = gpar(fontsize = text_size, fontface = "bold"),
                                          grid_width = unit(grid_width, "cm")),
                              Time=list(direction = "horizontal",
                                        nrow=1, 
                                        labels_gp = gpar(fontsize = text_size),
                                        title_gp = gpar(fontsize = text_size, fontface = "bold"),
                                        grid_width = unit(grid_width, "cm")))
                            ,col=list(Stim=colors,
                                      Time=dose_col),
                            annotation_name_gp = gpar(fontsize = text_size,fontface = "bold"),
                            simple_anno_size= unit(annot_height, "cm"))
col_annot@anno_list[["Time"]]@color_mapping@levels=names(dose_col)
col_annot@anno_list[["Time"]]@color_mapping@colors=dose_col

mean_scale= mean_sig[rownames(meta),]


row_ha = rowAnnotation(Coherence_score = anno_barplot(cs_list_ind[,"colli"], border = F,ylim = c(0,1),
                                                      axis_param = list(gp = gpar(fontsize = axis_size))  # Bold and larger font for axis labels
),
annotation_name_gp = gpar(fontsize = text_size, fontface = "bold"),
width = unit(2, "cm"))  # Bold and larger font for annotation name


colli_p=Heatmap(t(mean_scale),cluster_columns = T , right_annotation = row_ha,
                cluster_rows = T,name = "Mean signature score",
                top_annotation= col_annot,column_title = "Colli",
                show_column_names = F,
                column_title_gp  = gpar(fontsize = text_size, fontface = "bold"),
                row_names_side = "right",row_dend_side = "left",
                heatmap_legend_param = list(direction = "horizontal",
                                            legend_width = unit(20, "cm"),
                                            legend_gp  =gpar(fontsize = text_size),
                                            labels_gp  = gpar(fontsize = text_size),
                                            title_gp = gpar( fontsize = text_size)),
                row_names_max_width = max_text_width(
                  colnames(mean_scale), 
                  gp = gpar(fontsize=text_size)
                ),
                row_names_gp = gpar(col = c(rep("red",2), 
                                            rep("black",ncol(mean_scale)-2)),
                                    fontsize=text_size),
                row_dend_width = unit(3, "cm"),
                height = unit(40, "cm"),
                width = unit(20, "cm")
)


draw(colli_p, heatmap_legend_side = "top",
     annotation_legend_side = "bottom",newpage = T)

##### Rai ##### 
mean_sig= mean_sig_function(rai_c)

meta=rai_df
meta$time=meta$time
meta$time=as.numeric(meta$time)

meta=meta[order(meta$Stim),]


colors_select=colors[unique(meta$Stim)]

# dose_col=  c("dodgerblue4","dodgerblue","darkturquoise", "chartreuse", 
#              "chartreuse3","green",
#              "forestgreen")
dose_col=rev(brewer.pal(n=3, name="RdBu"))
names(dose_col)=sort(unique(meta$time))

colnames(meta)=gsub("time","Time",colnames(meta))

col_annot=HeatmapAnnotation(df = meta[,c("Stim","Time"),drop=F],
                            annotation_legend_param = list(
                              Stim = list(direction = "horizontal",
                                          nrow=1, 
                                          labels_gp = gpar(fontsize = text_size),
                                          title_gp = gpar(fontsize = text_size, fontface = "bold"),
                                          grid_width = unit(grid_width, "cm")),
                              Time=list(direction = "horizontal",
                                        nrow=1, 
                                        labels_gp = gpar(fontsize = text_size),
                                        title_gp = gpar(fontsize = text_size, fontface = "bold"),
                                        grid_width = unit(grid_width, "cm")))
                            ,col=list(Stim=colors,
                                      Time=dose_col),
                            annotation_name_gp = gpar(fontsize = text_size,fontface = "bold"),
                            simple_anno_size= unit(annot_height, "cm"))
col_annot@anno_list[["Time"]]@color_mapping@levels=names(dose_col)
col_annot@anno_list[["Time"]]@color_mapping@colors=dose_col

mean_scale= mean_sig[rownames(meta),]

row_ha = rowAnnotation(Coherence_score = anno_barplot(cs_list_ind[,"rai"], border = F,ylim = c(0,1),
                                                      axis_param = list(gp = gpar(fontsize = axis_size))  # Bold and larger font for axis labels
),
annotation_name_gp = gpar(fontsize = text_size, fontface = "bold"),
width = unit(2, "cm"))  # Bold and larger font for annotation name


rai_p=Heatmap(t(mean_scale),cluster_columns = T , right_annotation = row_ha,
              cluster_rows = T,name = "Mean signature score",
              top_annotation= col_annot,column_title = "Rai",
              show_column_names = F,
              column_title_gp  = gpar(fontsize = text_size, fontface = "bold"),
              row_names_side = "right",row_dend_side = "left",
              heatmap_legend_param = list(direction = "horizontal",
                                          legend_width = unit(20, "cm"),
                                          legend_gp  =gpar(fontsize = text_size),
                                          labels_gp  = gpar(fontsize = text_size),
                                          title_gp = gpar( fontsize = text_size)),
              row_names_max_width = max_text_width(
                colnames(mean_scale), 
                gp = gpar(fontsize=text_size)
              ),
              row_names_gp = gpar(col = c(rep("red",2), 
                                          rep("black",ncol(mean_scale)-2)),
                                  fontsize=text_size),
              width = unit(10, "cm"),
              row_dend_width = unit(3, "cm"),
              height = unit(40, "cm")
)

draw(rai_p, heatmap_legend_side = "top",
     annotation_legend_side = "bottom",newpage = T)

##### Fujiwara ##### 
mean_sig= mean_sig_function(fuji_c)

meta=fuji_df
meta=meta[order(meta$Stim),]

colors_select=colors[unique(meta$Stim)]

col_annot = HeatmapAnnotation(
  df = meta[, c("Stim"), drop = FALSE],
  annotation_legend_param = list(
    Stim = list(
      direction = "horizontal",
      nrow = 1,
      labels_gp = gpar(fontsize = text_size),
      title_gp = gpar(fontsize = text_size, fontface = "bold"),
      grid_width = unit(grid_width, "cm")
    )
  ),
  col = list(Stim = colors),
  annotation_name_gp = gpar(fontsize = text_size, fontface = "bold"),
  simple_anno_size= unit(annot_height, "cm")
)


mean_scale= mean_sig[rownames(meta),]


row_ha = rowAnnotation(Coherence_score = anno_barplot(cs_list_ind[,"fujiwara"], border = F,ylim = c(0,1),
                                                      axis_param = list(gp = gpar(fontsize = axis_size))  # Bold and larger font for axis labels
),
annotation_name_gp = gpar(fontsize = text_size, fontface = "bold"),
width = unit(2, "cm"))  # Bold and larger font for annotation name


fuji_p=Heatmap(t(mean_scale),cluster_columns = T , right_annotation = row_ha,
               cluster_rows = T,name = "Mean signature score",
               top_annotation= col_annot,column_title = "Fujiwara",
               column_title_gp  = gpar(fontsize = text_size, fontface = "bold"),
               show_column_names = F,
               row_names_side = "right",row_dend_side = "left",
               heatmap_legend_param = list(direction = "horizontal",
                                           legend_width = unit(20, "cm"),
                                           legend_gp  =gpar(fontsize = text_size),
                                           labels_gp  = gpar(fontsize = text_size),
                                           title_gp = gpar( fontsize = text_size),
                                           position = "top"),
               row_names_max_width = max_text_width(
                 colnames(mean_scale), 
                 gp = gpar(fontsize=text_size)
               ),
               row_names_gp = gpar(col = c(rep("red",2), 
                                           rep("black",ncol(mean_scale)-2)),
                                   fontsize=text_size),
               width = unit(10, "cm"),
               height = unit(40, "cm"),
               row_dend_width = unit(3, "cm")
)

draw(fuji_p, heatmap_legend_side = "top",
     annotation_legend_side = "bottom",newpage = T)

fuji_grob <- grid::grid.grabExpr(draw(fuji_p, heatmap_legend_side = "top",
                                         annotation_legend_side = "top",newpage = T))
rai_grob <- grid::grid.grabExpr(draw(rai_p, heatmap_legend_side = "top",
                                         annotation_legend_side = "top",newpage = T))
jank_grob <- grid::grid.grabExpr(draw(jank_p, heatmap_legend_side = "top",
                                      annotation_legend_side = "top",newpage = T))
colli_grob <- grid::grid.grabExpr(draw(colli_p, heatmap_legend_side = "top",
                                         annotation_legend_side = "top",newpage = T))
ziegler_grob <- grid::grid.grabExpr(draw(ziegler_beas_p, heatmap_legend_side = "top",
                                         annotation_legend_side = "top",newpage = T))



# Arrange the heatmap and ggplot together

pdf(paste0(save_dir,"Discovery_heatmap_all.pdf"),
    width = 45,height = 50 )
plot_grid(
  plot_grid(jank_grob,
            fuji_grob,rai_grob, nrow = 1),
  plot_grid(ziegler_grob,colli_grob,ncol=2),nrow=2,
  labels = c("A",""),vjust = 2,label_size = 100)


dev.off()
