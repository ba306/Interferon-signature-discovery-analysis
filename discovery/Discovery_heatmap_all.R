# Heatmaps for our IFN-I and IFN-II signatures along with published IFN signatures
# In discovery datasets

# Load required libraries
library(foreach)
library(dplyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(gridtext)
library(grid)
library(cowplot)
source("functions/CScores.R") # coherence score calculation

# Predefined parameters
text_size    = 16   # Increased for better visibility
axis_size    = 1
grid_width   = 0.5
annot_height = 0.5
heat_width   = 12   # Increased width for better spacing
heat_height  = 14
row_denth    = 1
legend_width = 5
CS_width     = 1
legend_height= 1.5  # Adjusted height for the color scale (legend)

# Directories
save_dir   = "discovery/" 
data_dir   = "data/IFN_bulk/" 

### Load IFN gene signatures ###
rosetta_new= readRDS(file = "published_signatures/IFN_signature_list_all_valid.RDS")


## Colors for stimulations 
colors = sapply(c("#e41a1c",  # red
                  "#377eb8",  # blue
                  "#4daf4a",  # green
                  "#984ea3",  # purple
                  "#ff7f00",  # orange
                  "#ffff33",  # yellow
                  "#a65628"), toupper)
colors = c(colors[1], colors[2], colors[3], colors[4], colors[5])
names(colors) = c("IFNa", "IFNb", "IFNg", "Control", "IFNlambda")


### Load discovery expression and metadata ###
load(file = paste0(data_dir, "discovery_exp_meta_rdata.RData"))

### Coherence score calculation ###
# For fast calculation, only genes from our gene lists are used
gene_coherence = unlist(rosetta_new, use.names = FALSE)

exp_list = list(ziegler_c, jank_c, colli_c, rai_c, fuji_c)

cs_list_ind = foreach(j = seq_along(exp_list), .combine = cbind) %do% {
  print(j)
  foreach(i = seq_along(rosetta_new), .combine = rbind) %do% {
    print(i)
    gen = intersect(rownames(exp_list[[j]]), rosetta_new[[i]])
    pear = cor(t(exp_list[[j]][gen, ]), method = "pearson")
    score = coherence_score_val(gen, pear)
    round(score, digits = 2)
  }
}
rownames(cs_list_ind) = names(rosetta_new)
colnames(cs_list_ind) = c("ziegler", "jankowski", "colli", "rai", "fujiwara")
# write.table(cs_list_ind, paste0(save_dir,"valid_cscores.tsv"), sep = "\t")

### Mean signature function ###
mean_sig_function = function(data) {
  data = as.matrix(data)
  data = t(scale(t(data), center = TRUE, scale = TRUE))
  
  foreach(i = seq_along(names(rosetta_new)), .combine = cbind) %do% {
    print(i)
    inter_genes = intersect(rosetta_new[[i]], rownames(data))
    data[inter_genes, , drop = FALSE] %>%
      colMeans(na.rm = TRUE) %>%
      as.data.frame() %>%
      setNames(nm = names(rosetta_new)[i])
  }
}

# Replace "untreated" with "Control" in metadata for all relevant data frames
var_names = names(as.list(.GlobalEnv))
df_names  = var_names[grepl("(.)_df", var_names)]
for(i in df_names){
  var = get(i)
  var$Stim = gsub("untreated", "Control", var$Stim)
  assign(x = i, var)
}

##########################
### ZIEGLER DONOR 1  #####
##########################
mean_sig = mean_sig_function(ziegler_c)
meta = ziegler_df
meta$Dose = as.numeric(meta$Dose)
# Order meta data by Stim and Dose (ascending)
meta = meta[order(meta$Stim, meta$Dose), ]

colors_select = colors[unique(meta$Stim)]
dose_col = rev(brewer.pal(n = 7, name = "RdBu"))
names(dose_col) = sort(unique(meta$Dose))

col_annot = HeatmapAnnotation(
  df = meta[, c("Stim", "Dose"), drop = FALSE],
  annotation_legend_param = list(
    Stim = list(
      direction = "vertical",
      nrow = 1, 
      labels_gp = gpar(fontsize = text_size),
      title_gp  = gpar(fontsize = text_size, fontface = "bold"),
      grid_width = unit(grid_width, "cm")
    ),
    Dose = list(
      direction = "vertical",
      nrow = 3, 
      labels_gp = gpar(fontsize = text_size),
      title_gp  = gpar(fontsize = text_size, fontface = "bold"),
      grid_width = unit(grid_width, "cm")
    )
  ),
  col = list(
    Stim = colors,
    Dose = dose_col
  ),
  annotation_name_gp = gpar(fontsize = text_size, fontface = "bold"),
  simple_anno_size = unit(annot_height, "cm")
)
col_annot@anno_list[["Dose"]]@color_mapping@levels = names(dose_col)
col_annot@anno_list[["Dose"]]@color_mapping@colors  = dose_col

mean_scale = mean_sig[rownames(meta), ]

row_ha = rowAnnotation(
  CS = anno_barplot(
    cs_list_ind[, "ziegler"],
    border = FALSE,
    ylim = c(0, 1),
    axis_param = list(
      gp = gpar(fontsize = 14),        # Increased fontsize for y-axis labels
      at = c(0, 0.5, 1),
      labels = c("0", "0.5", "1")
    ),
    labels = cs_list_ind[, "ziegler"],
    label_gp = gpar(fontsize = 100)    # Increased fontsize for bar labels
  ),
  annotation_name_gp = gpar(fontsize = text_size, fontface = "bold"),
  width = unit(CS_width, "cm")
)

ziegler_beas_p = Heatmap(
  t(mean_scale),
  cluster_columns = TRUE,
  right_annotation = row_ha,
  cluster_rows = TRUE,
  name = "Mean signature score",
  top_annotation = col_annot,
  column_title = "Ziegler Donor 1",
  column_title_gp = gpar(fontsize = text_size, fontface = "bold"),
  show_column_names = FALSE,
  row_names_side = "right",
  row_dend_side = "left",
  heatmap_legend_param = list(
    direction = "horizontal",
    legend_width = unit(legend_width, "cm"),
    legend_height = unit(legend_height, "cm"),
    legend_gp = gpar(fontsize = text_size),
    labels_gp = gpar(fontsize = text_size),
    title_gp = gpar(fontsize = text_size)
  ),
  row_names_max_width = max_text_width(colnames(mean_scale), gp = gpar(fontsize = text_size)),
  row_names_gp = gpar(
    col = c(rep("red", 2), rep("black", ncol(mean_scale) - 2)),
    fontsize = text_size
  ),
  row_dend_width = unit(row_denth, "cm"),
  height = unit(heat_height, "cm"),
  width = unit(heat_width, "cm")
)

draw(ziegler_beas_p, heatmap_legend_side = "top", annotation_legend_side = "bottom", newpage = TRUE)
ziegler_grob = grid::grid.grabExpr(draw(ziegler_beas_p,
                                        heatmap_legend_side = "top",
                                        annotation_legend_side = "bottom",
                                        newpage = TRUE))

##########################
### JANKOWSKI         ###
##########################
mean_sig = mean_sig_function(jank_c)
meta = jank_df
# Order meta if needed (by Stim)
meta = meta[order(meta$Stim), ]
colors_select = colors[unique(meta$Stim)]

col_annot = HeatmapAnnotation(
  df = meta[, c("Stim"), drop = FALSE],
  annotation_legend_param = list(
    Stim = list(
      direction = "horizontal",
      nrow = 2,
      grid_width = unit(grid_width, "cm"),
      labels_gp = gpar(fontsize = text_size),
      title_gp  = gpar(fontsize = text_size, fontface = "bold")
    )
  ),
  col = list(Stim = colors),
  annotation_name_gp = gpar(fontsize = text_size, fontface = "bold"),
  simple_anno_size = unit(annot_height, "cm")
)
mean_scale = mean_sig[rownames(meta), ]

row_ha = rowAnnotation(
  CS = anno_barplot(
    cs_list_ind[, "jankowski"],
    border = FALSE,
    ylim = c(0, 1),
    axis_param = list(
      gp = gpar(fontsize = 14),
      at = c(0, 0.5, 1),
      labels = c("0", "0.5", "1")
    )
  ),
  annotation_name_gp = gpar(fontsize = text_size, fontface = "bold"),
  width = unit(CS_width, "cm")
)

jank_p = Heatmap(
  t(mean_scale),
  cluster_columns = TRUE,
  right_annotation = row_ha,
  cluster_rows = TRUE,
  name = "Mean signature score",
  top_annotation = col_annot,
  column_title = "Jankowski",
  column_title_gp = gpar(fontsize = text_size, fontface = "bold"),
  show_column_names = FALSE,
  row_names_side = "right",
  row_dend_side = "left",
  heatmap_legend_param = list(
    direction = "horizontal",
    legend_width = unit(legend_width, "cm"),
    legend_gp = gpar(fontsize = text_size),
    labels_gp = gpar(fontsize = text_size),
    title_gp = gpar(fontsize = text_size)
  ),
  row_names_max_width = max_text_width(colnames(mean_scale), gp = gpar(fontsize = text_size)),
  row_names_gp = gpar(
    col = c(rep("red", 2), rep("black", ncol(mean_scale) - 2)),
    fontsize = text_size
  ),
  width = unit(heat_width/2, "cm"),
  row_dend_width = unit(row_denth, "cm"),
  height = unit(heat_height, "cm")
)

draw(jank_p, heatmap_legend_side = "top", annotation_legend_side = "bottom", newpage = TRUE)
jank_grob = grid::grid.grabExpr(draw(jank_p,
                                     heatmap_legend_side = "top",
                                     annotation_legend_side = "bottom",
                                     newpage = TRUE))

##########################
### COLLI              ###
##########################
mean_sig = mean_sig_function(colli_c)
meta = colli_df
meta$time = as.numeric(meta$time)
# Order meta data by Stim and time (ascending)
meta = meta[order(meta$Stim, meta$time), ]
colors_select = colors[unique(meta$Stim)]
dose_col = rev(brewer.pal(n = 3, name = "RdBu"))
names(dose_col) = sort(unique(meta$time))
colnames(meta) = gsub("time", "Time", colnames(meta))

col_annot = HeatmapAnnotation(
  df = meta[, c("Stim", "Time"), drop = FALSE],
  annotation_legend_param = list(
    Stim = list(
      direction = "horizontal",
      nrow = 1,
      labels_gp = gpar(fontsize = text_size),
      title_gp  = gpar(fontsize = text_size, fontface = "bold"),
      grid_width = unit(grid_width, "cm")
    ),
    Time = list(
      direction = "horizontal",
      nrow = 1,
      labels_gp = gpar(fontsize = text_size),
      title_gp  = gpar(fontsize = text_size, fontface = "bold"),
      grid_width = unit(grid_width, "cm")
    )
  ),
  col = list(
    Stim = colors,
    Time = dose_col
  ),
  annotation_name_gp = gpar(fontsize = text_size, fontface = "bold"),
  simple_anno_size = unit(annot_height, "cm")
)
col_annot@anno_list[["Time"]]@color_mapping@levels = names(dose_col)
col_annot@anno_list[["Time"]]@color_mapping@colors  = dose_col

mean_scale = mean_sig[rownames(meta), ]

row_ha = rowAnnotation(
  CS = anno_barplot(
    cs_list_ind[, "colli"],
    border = FALSE,
    ylim = c(0, 1),
    axis_param = list(
      gp = gpar(fontsize = 14),
      at = c(0, 0.5, 1),
      labels = c("0", "0.5", "1")
    )
  ),
  annotation_name_gp = gpar(fontsize = text_size, fontface = "bold"),
  width = unit(CS_width, "cm")
)

colli_p = Heatmap(
  t(mean_scale),
  cluster_columns = TRUE,
  right_annotation = row_ha,
  cluster_rows = TRUE,
  name = "Mean signature score",
  top_annotation = col_annot,
  column_title = "Colli",
  column_title_gp = gpar(fontsize = text_size, fontface = "bold"),
  show_column_names = FALSE,
  row_names_side = "right",
  row_dend_side = "left",
  heatmap_legend_param = list(
    direction = "horizontal",
    legend_width = unit(legend_width, "cm"),
    legend_gp = gpar(fontsize = text_size),
    labels_gp = gpar(fontsize = text_size),
    title_gp = gpar(fontsize = text_size)
  ),
  row_names_max_width = max_text_width(colnames(mean_scale), gp = gpar(fontsize = text_size)),
  row_names_gp = gpar(
    col = c(rep("red", 2), rep("black", ncol(mean_scale) - 2)),
    fontsize = text_size
  ),
  row_dend_width = unit(row_denth, "cm"),
  height = unit(heat_height, "cm"),
  width = unit(heat_width, "cm")
)

draw(colli_p, heatmap_legend_side = "top", annotation_legend_side = "bottom", newpage = TRUE)
colli_grob = grid::grid.grabExpr(draw(colli_p,
                                      heatmap_legend_side = "top",
                                      annotation_legend_side = "bottom",
                                      newpage = TRUE))

##########################
### RAI                ###
##########################
mean_sig = mean_sig_function(rai_c)
meta = rai_df
meta$time = as.numeric(meta$time)
# Order meta data by Stim and time (ascending)
meta = meta[order(meta$Stim, meta$time), ]
colors_select = colors[unique(meta$Stim)]
dose_col = rev(brewer.pal(n = 3, name = "RdBu"))
names(dose_col) = sort(unique(meta$time))
colnames(meta) = gsub("time", "Time", colnames(meta))

col_annot = HeatmapAnnotation(
  df = meta[, c("Stim", "Time"), drop = FALSE],
  annotation_legend_param = list(
    Stim = list(
      direction = "vertical",
      ncol = 1,
      labels_gp = gpar(fontsize = text_size),
      title_gp  = gpar(fontsize = text_size, fontface = "bold"),
      grid_width = unit(grid_width, "cm")
    ),
    Time = list(
      direction = "vertical",
      ncol = 1,
      labels_gp = gpar(fontsize = text_size),
      title_gp  = gpar(fontsize = text_size, fontface = "bold"),
      grid_width = unit(grid_width, "cm")
    )
  ),
  col = list(
    Stim = colors,
    Time = dose_col
  ),
  annotation_name_gp = gpar(fontsize = text_size, fontface = "bold"),
  simple_anno_size = unit(annot_height, "cm")
)
col_annot@anno_list[["Time"]]@color_mapping@levels = names(dose_col)
col_annot@anno_list[["Time"]]@color_mapping@colors  = dose_col

mean_scale = mean_sig[rownames(meta), ]

row_ha = rowAnnotation(
  CS = anno_barplot(
    cs_list_ind[, "rai"],
    border = FALSE,
    ylim = c(0, 1),
    axis_param = list(
      gp = gpar(fontsize = 14),
      at = c(0, 0.5, 1),
      labels = c("0", "0.5", "1")
    )
  ),
  annotation_name_gp = gpar(fontsize = text_size, fontface = "bold"),
  width = unit(CS_width, "cm")
)

rai_p = Heatmap(
  t(mean_scale),
  cluster_columns = TRUE,
  right_annotation = row_ha,
  cluster_rows = TRUE,
  name = "Mean signature score",
  top_annotation = col_annot,
  column_title = "Rai",
  column_title_gp = gpar(fontsize = text_size, fontface = "bold"),
  show_column_names = FALSE,
  row_names_side = "right",
  row_dend_side = "left",
  heatmap_legend_param = list(
    direction = "horizontal",
    legend_width = unit(legend_width, "cm"),
    legend_gp = gpar(fontsize = text_size),
    labels_gp = gpar(fontsize = text_size),
    title_gp = gpar(fontsize = text_size)
  ),
  row_names_max_width = max_text_width(colnames(mean_scale), gp = gpar(fontsize = text_size)),
  row_names_gp = gpar(
    col = c(rep("red", 2), rep("black", ncol(mean_scale) - 2)),
    fontsize = text_size
  ),
  width = unit(heat_width/2, "cm"),
  row_dend_width = unit(row_denth, "cm"),
  height = unit(heat_height, "cm")
)

draw(rai_p, heatmap_legend_side = "top", annotation_legend_side = "bottom", newpage = TRUE)
rai_grob = grid::grid.grabExpr(draw(rai_p,
                                    heatmap_legend_side = "top",
                                    annotation_legend_side = "bottom",
                                    newpage = TRUE))

##########################
### FUJIWARA           ###
##########################
mean_sig = mean_sig_function(fuji_c)
meta = fuji_df
# Order meta if needed (by Stim)
meta = meta[order(meta$Stim), ]
colors_select = colors[unique(meta$Stim)]

col_annot = HeatmapAnnotation(
  df = meta[, c("Stim"), drop = FALSE],
  annotation_legend_param = list(
    Stim = list(
      direction = "horizontal",
      nrow = 1,
      labels_gp = gpar(fontsize = text_size),
      title_gp  = gpar(fontsize = text_size, fontface = "bold"),
      grid_width = unit(grid_width, "cm")
    )
  ),
  col = list(Stim = colors),
  annotation_name_gp = gpar(fontsize = text_size, fontface = "bold"),
  simple_anno_size = unit(annot_height, "cm")
)
mean_scale = mean_sig[rownames(meta), ]

row_ha = rowAnnotation(
  CS = anno_barplot(
    cs_list_ind[, "fujiwara"],
    border = FALSE,
    ylim = c(0, 1),
    axis_param = list(
      gp = gpar(fontsize = 14),
      at = c(0, 0.5, 1),
      labels = c("0", "0.5", "1")
    )
  ),
  annotation_name_gp = gpar(fontsize = text_size, fontface = "bold"),
  width = unit(CS_width, "cm")
)

fuji_p = Heatmap(
  t(mean_scale),
  cluster_columns = TRUE,
  right_annotation = row_ha,
  cluster_rows = TRUE,
  name = "Mean signature score",
  top_annotation = col_annot,
  column_title = "Fujiwara",
  column_title_gp = gpar(fontsize = text_size, fontface = "bold"),
  show_column_names = FALSE,
  row_names_side = "right",
  row_dend_side = "left",
  heatmap_legend_param = list(
    direction = "horizontal",
    legend_width = unit(legend_width, "cm"),
    legend_gp = gpar(fontsize = text_size),
    labels_gp = gpar(fontsize = text_size),
    title_gp = gpar(fontsize = text_size),
    position = "top"
  ),
  row_names_max_width = max_text_width(colnames(mean_scale), gp = gpar(fontsize = text_size)),
  row_names_gp = gpar(
    col = c(rep("red", 2), rep("black", ncol(mean_scale) - 2)),
    fontsize = text_size
  ),
  width = unit(heat_width/2, "cm"),
  height = unit(heat_height, "cm"),
  row_dend_width = unit(row_denth, "cm")
)

draw(fuji_p, heatmap_legend_side = "top", annotation_legend_side = "bottom", newpage = TRUE)
fuji_grob = grid::grid.grabExpr(draw(fuji_p,
                                     heatmap_legend_side = "top",
                                     annotation_legend_side = "bottom",
                                     newpage = TRUE))

##########################
### Arrangement of all ###
##########################
pdf(paste0(save_dir,"Discovery_heatmap_all_res.pdf"),
    width = 18,height = 18 )
plot_grid(
  plot_grid(jank_grob,
            fuji_grob,rai_grob, nrow = 1),
  plot_grid(ziegler_grob,colli_grob,ncol=2),nrow=2,
  labels = c("A",""),vjust = 2,label_size = 26)


dev.off()

