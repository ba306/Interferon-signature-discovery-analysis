# Heatmaps for our IFN-I and IFN-II signatures along with published IFN signatures
# In validation datasets

# Load required libraries
library(foreach)
library(dplyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(gridtext)
library(grid)
library(cowplot) 
source("functions/CScores.R") # calculating Coherence score

# Predefined heatmap parameters
text_size    = 16   # Increased for better visibility
axis_size    = 1
grid_width   = 0.5
annot_height = 0.5
heat_width   = 12   # Width for heatmap (in cm)
heat_height  = 14   # Height for heatmap (in cm)
row_denth    = 1    # Row dendrogram width (in cm)
legend_width = 5
CS_width     = 1
legend_height= 1.5  # Legend height (in cm)

# Directories
save_dir   = "discovery/" 
data_dir   = "data/IFN_bulk/" 


# Load validation expression and metadata  

load(
  file = paste0(data_dir,
                "validation_exp_meta_rdata.RData"))

### Load IFN gene signatures ###
rosetta_new = readRDS(file = "published_signatures/IFN_signature_list_all_valid.RDS")

## Colors for stimulations 
colors = sapply(c("#e41a1c",  "#377eb8",  "#4daf4a",  "#984ea3",  "#ff7f00",  "#ffff33",  "#a65628"), toupper)
colors = c(colors[1], colors[2], colors[3], colors[4], colors[5])
names(colors) = c("IFNa", "IFNb", "IFNg", "Control", "IFNlambda")

### Coherence score calculation ###

# For fast calculation only genes from our gene lists used
gene_coherence = unlist(rosetta_new, use.names = FALSE)
exp_list = list(ziegler_donor2_list[[1]],
                ziegler_beas_list[[1]],
                devlin_list[[1]],
                lee_list[[1]])

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
colnames(cs_list_ind) = c("ziegler_donor2", "ziegler_BEAS2B", "devlin", "lee")

### Individual heatmaps ####

# Function to calculate mean signature
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

#### Ziegler BEAS-2B ####
mean_sig = mean_sig_function(ziegler_beas_list[[1]])
meta = ziegler_beas_list[[2]]
meta$Stim = gsub("IFNA", "IFNa", meta$Stim)
meta$Stim = gsub("IFNG", "IFNg", meta$Stim)
meta$Stim = gsub("Untreated", "Control", meta$Stim)
meta$Dose = as.numeric(meta$Dose)
meta = meta[order(meta$Stim, meta$Dose), ]

colors_select = colors[unique(meta$Stim)]
dose_col = rev(brewer.pal(n = 7, name = "RdBu"))
names(dose_col) = sort(unique(meta$Dose))

col_annot = HeatmapAnnotation(
  df = meta[, c("Stim", "Dose"), drop = FALSE],
  annotation_legend_param = list(
    Stim = list(direction = "vertical",
                nrow = 1, 
                labels_gp = gpar(fontsize = text_size),
                title_gp = gpar(fontsize = text_size, fontface = "bold"),
                grid_width = unit(grid_width, "cm")),
    Dose = list(direction = "vertical",
                nrow = 3, 
                labels_gp = gpar(fontsize = text_size),
                title_gp = gpar(fontsize = text_size, fontface = "bold"),
                grid_width = unit(grid_width, "cm"))
  ),
  col = list(Stim = colors, Dose = dose_col),
  annotation_name_gp = gpar(fontsize = text_size, fontface = "bold"),
  simple_anno_size = unit(annot_height, "cm")
)
col_annot@anno_list[["Dose"]]@color_mapping@levels = names(dose_col)
col_annot@anno_list[["Dose"]]@color_mapping@colors = dose_col

mean_scale = mean_sig[rownames(meta), ]

row_ha = rowAnnotation(
  CS = anno_barplot(
    cs_list_ind[, "ziegler_BEAS2B"],
    border = FALSE,
    ylim = c(0, 1),
    axis_param = list(
      gp = gpar(fontsize = 14),        # Increased fontsize for y-axis labels
      at = c(0, 0.5, 1),
      labels = c("0", "0.5", "1")
    ),
    labels = cs_list_ind[, "ziegler_BEAS2B"],
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
  column_title = "Ziegler BEAS-2B",
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

#### Devlin ####
mean_sig = mean_sig_function(devlin_list[[1]])
meta = devlin_list[[2]]
meta$Stim = gsub("IFNbeta", "IFNb", meta$Stim)
meta$Stim = gsub("IFNgamma", "IFNg", meta$Stim)
meta$Stim = gsub("Buffer", "Control", meta$Stim)
meta = meta[order(meta$Stim), ]

colors_select = colors[unique(meta$Stim)]
col_annot = HeatmapAnnotation(
  df = meta[, c("Stim"), drop = FALSE],
  annotation_legend_param = list(
    Stim = list(direction = "horizontal",
                grid_width = unit(grid_width, "cm"),
                nrow = 1, 
                labels_gp = gpar(fontsize = text_size),
                title_gp = gpar(fontsize = text_size, fontface = "bold"))
  ),
  col = list(Stim = colors),
  annotation_name_gp = gpar(fontsize = text_size, fontface = "bold"),
  simple_anno_size = unit(annot_height, "cm")
)
mean_scale = mean_sig[rownames(meta), ]
row_ha = rowAnnotation(
  CS = anno_barplot(
    cs_list_ind[, "devlin"],
    border = FALSE,
    ylim = c(0, 1),
    axis_param = list(gp = gpar(fontsize = 14),
                      at = c(0, 0.5, 1),
                      labels = c("0", "0.5", "1")),
    labels = cs_list_ind[, "devlin"],
    label_gp = gpar(fontsize = 100)
  ),
  annotation_name_gp = gpar(fontsize = text_size, fontface = "bold"),
  width = unit(CS_width, "cm")
)
devlin_p = Heatmap(
  t(mean_scale), 
  cluster_columns = TRUE, 
  right_annotation = row_ha,
  cluster_rows = TRUE,
  name = "Mean signature score",
  top_annotation = col_annot,
  column_title = "Devlin",
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
  width = unit(heat_width, "cm"),
  row_dend_width = unit(row_denth, "cm"),
  height = unit(heat_height, "cm")
)

draw(devlin_p, heatmap_legend_side = "top", annotation_legend_side = "bottom", newpage = TRUE)

#### Ziegler Donor 2 ####
mean_sig = mean_sig_function(ziegler_donor2_list[[1]])
meta = ziegler_donor2_list[[2]]
meta$Stim = gsub("IFNA", "IFNa", meta$Stim)
meta$Stim = gsub("IFNG", "IFNg", meta$Stim)
meta$Stim = gsub("Untreated", "Control", meta$Stim)
meta$Dose = as.numeric(meta$Dose)
meta = meta[order(meta$Stim), ]

colors_select = colors[unique(meta$Stim)]
dose_col = rev(brewer.pal(n = 7, name = "RdBu"))
names(dose_col) = sort(unique(meta$Dose))
col_annot = HeatmapAnnotation(
  df = meta[, c("Stim", "Dose"), drop = FALSE],
  annotation_legend_param = list(
    Stim = list(direction = "horizontal",
                nrow = 1, 
                labels_gp = gpar(fontsize = text_size),
                title_gp = gpar(fontsize = text_size, fontface = "bold"),
                grid_width = unit(grid_width, "cm")),
    Dose = list(direction = "horizontal",
                nrow = 3, 
                labels_gp = gpar(fontsize = text_size),
                title_gp = gpar(fontsize = text_size, fontface = "bold"),
                grid_width = unit(grid_width, "cm"))
  ),
  col = list(Stim = colors, Dose = dose_col),
  annotation_name_gp = gpar(fontsize = text_size, fontface = "bold"),
  simple_anno_size = unit(annot_height, "cm")
)
col_annot@anno_list[["Dose"]]@color_mapping@levels = names(dose_col)
col_annot@anno_list[["Dose"]]@color_mapping@colors = dose_col
mean_scale = mean_sig[rownames(meta), ]
row_ha = rowAnnotation(
  CS = anno_barplot(
    cs_list_ind[, "ziegler_donor2"],
    border = FALSE,
    ylim = c(0, 1),
    axis_param = list(gp = gpar(fontsize = 14),
                      at = c(0, 0.5, 1),
                      labels = c("0", "0.5", "1")),
    labels = cs_list_ind[, "ziegler_donor2"],
    label_gp = gpar(fontsize = 100)
  ),
  annotation_name_gp = gpar(fontsize = text_size, fontface = "bold"),
  width = unit(CS_width, "cm")
)
ziegler_2_p = Heatmap(
  t(mean_scale), 
  cluster_columns = TRUE, 
  right_annotation = row_ha,
  cluster_rows = TRUE,
  name = "Mean signature score",
  top_annotation = col_annot,
  column_title = "Ziegler Donor 2",
  show_column_names = FALSE,
  column_title_gp = gpar(fontsize = text_size, fontface = "bold"),
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

draw(ziegler_2_p, heatmap_legend_side = "top", annotation_legend_side = "bottom", newpage = TRUE)

#### Lee ####
mean_sig = mean_sig_function(lee_list[[1]])
meta = lee_list[[2]]
meta$Stim = gsub("IFNl3", "IFNlambda", meta$Stim)
meta = meta[order(meta$Stim), ]
colors_select = colors[unique(meta$Stim)]
col_annot = HeatmapAnnotation(
  df = meta[, c("Stim"), drop = FALSE],
  annotation_legend_param = list(
    Stim = list(direction = "horizontal",
                grid_width = unit(grid_width, "cm"),
                nrow = 1, 
                labels_gp = gpar(fontsize = text_size),
                title_gp = gpar(fontsize = text_size, fontface = "bold"))
  ),
  col = list(Stim = colors),
  annotation_name_gp = gpar(fontsize = text_size, fontface = "bold"),
  simple_anno_size = unit(annot_height, "cm")
)
mean_scale = mean_sig[rownames(meta), ]
row_ha = rowAnnotation(
  CS = anno_barplot(
    cs_list_ind[, "lee"],
    border = FALSE,
    ylim = c(0, 1),
    axis_param = list(gp = gpar(fontsize = 14),
                      at = c(0, 0.5, 1),
                      labels = c("0", "0.5", "1")),
    labels = cs_list_ind[, "lee"],
    label_gp = gpar(fontsize = 100)
  ),
  annotation_name_gp = gpar(fontsize = text_size, fontface = "bold"),
  width = unit(CS_width, "cm")
)
lee_p = Heatmap(
  t(mean_scale), 
  cluster_columns = TRUE, 
  right_annotation = row_ha,
  cluster_rows = TRUE,
  name = "Mean signature score",
  top_annotation = col_annot,
  column_title = "Lee",
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
  width = unit(heat_width, "cm"),
  row_dend_width = unit(row_denth, "cm"),
  height = unit(heat_height, "cm")
)

draw(lee_p, heatmap_legend_side = "top", annotation_legend_side = "bottom", newpage = TRUE)

#### Putting all plots together ####
ziegler_2_p_grob <- grid::grid.grabExpr(draw(ziegler_2_p, heatmap_legend_side = "top",
                                             annotation_legend_side = "bottom", newpage = TRUE))
ziegler_beas_p_grob <- grid::grid.grabExpr(draw(ziegler_beas_p, heatmap_legend_side = "top",
                                                annotation_legend_side = "bottom", newpage = TRUE))
devlin_p_grob <- grid::grid.grabExpr(draw(devlin_p, heatmap_legend_side = "top",
                                          annotation_legend_side = "bottom", newpage = TRUE))
lee_p_grob <- grid::grid.grabExpr(draw(lee_p, heatmap_legend_side = "top",
                                       annotation_legend_side = "bottom", newpage = TRUE))

pdf(paste0(save_dir, "Validation_heatmap_all_res.pdf"),  width = 17,height = 17)
plot_grid(
  ziegler_2_p_grob, ziegler_beas_p_grob, devlin_p_grob, lee_p_grob, nrow = 2,
  labels = c("B", "", "", ""), vjust = 2, label_size = 26
)
dev.off()

