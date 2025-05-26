# Function to create heatmap for a given signature
# Choose top correlated or covariated signatures

rosettasx_heatmap_func=function(data,
                                covariance=F,
                                correlation=F,
                                selected_signature=NULL,
                                top=20, # top 20 correlated or covariance 
                                textsize=8, # size for top signature labels
                                # titlesize=8, 
                                label_size=12, # size for selected signature labels
                                signature_colors=NULL,# signature score colors
                                row_annot_colors=NULL # covariance or correlation color scale
){
  
  
  # order samples by selected signature
  col_order=order(data[,selected_signature])
  
  # Create column annotations for main heatmap
  exp_selected=data[,selected_signature,drop=F]
  colnames(exp_selected)="exp"
  
  col_anno=ComplexHeatmap::HeatmapAnnotation(
    df                   = exp_selected,
    show_legend          = F,
    na_col               = "grey",
    show_annotation_name = TRUE,
    col=list(exp=signature_colors),
    simple_anno_size_adjust  = TRUE,
    annotation_label = selected_signature)
  
  
  # Either correlation or covariance based column annotation
  
  if(covariance==T){
    # Covariance between selected signature and other signatures 
    cov_Df=data %>%
      cov() %>%
      .[, grepl(selected_signature, colnames(data))] %>%
      as.matrix()
    
    cov_Df=cov_Df[!grepl(selected_signature,rownames(cov_Df)),,drop=F]
    
    colnames(cov_Df)="Covariance"
    
    top_signatures=order(cov_Df[,1],decreasing = T)[1:top]
    
    row_anno=ComplexHeatmap::HeatmapAnnotation(
      df                   = cov_Df[top_signatures,,drop=F],
      show_legend          = T,
      na_col               = "grey",
      show_annotation_name = TRUE,
      col=list(Covariance=row_annot_colors),
      simple_anno_size_adjust  = TRUE,
      which="row",
      annotation_legend_param  = base::list(
        direction      = "horizontal",
        # title_position = "topleft",
        title_gp = gpar(fontsize = label_size, fontface = "bold"),
        labels_gp = gpar(fontsize = label_size)
        
      )
    )
    
  }else if(correlation==T){
    # Correlation between selected signature and other signatures 
    cor_Df=data %>%
      cor() %>%
      .[, grepl(selected_signature, colnames(data))] %>%
      as.matrix()
    
    cor_Df=cor_Df[!grepl(selected_signature,rownames(cor_Df)),,drop=F]
    
    colnames(cor_Df)="Correlation"
    
    top_signatures=order(cor_Df[,1],decreasing = T)[1:top]
    
    row_anno=ComplexHeatmap::HeatmapAnnotation(
      df                   = cor_Df[top_signatures,,drop=F],
      show_legend          = T,
      na_col               = "grey",
      show_annotation_name = TRUE,
      col=list(Correlation=row_annot_colors),
      simple_anno_size_adjust  = TRUE,
      which="row",
      annotation_legend_param  = base::list(
        direction      = "horizontal",
        # title_position = "topleft",
        title_gp = gpar(fontsize = label_size, fontface = "bold"),
        labels_gp = gpar(fontsize = label_size)
        
      )
    )
    
  }
  
  # Plotting
  # Main heatmap
  
  data_heatmap=data[,!grepl(selected_signature,colnames(data))]
  
  main_ht=Heatmap(t(data_heatmap[,  top_signatures]),show_column_names = F,cluster_rows = F,
                  cluster_columns = F,
                  column_order = col_order,
                  top_annotation       = col_anno,
                  left_annotation  = row_anno,
                  col                  = signature_colors,
                  heatmap_legend_param = base::list(
                    direction      = "horizontal",
                    # title_position = "topleft",
                    title_gp = gpar(fontsize = label_size, fontface = "bold"),
                    labels_gp = gpar(fontsize = label_size)
                    
                  ),name="Signature score",row_names_gp = gpar(fontsize = textsize)
                  
  )
  draw(main_ht,heatmap_legend_side="top",
       annotation_legend_side      = "bottom"
  )
  
}
