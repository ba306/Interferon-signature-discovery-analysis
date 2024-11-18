# Prediction metrics calculation

library(varhandle)
library(foreach)
library(dplyr)
set.seed(42)


##### Calculate prediction metrics ##### 

# Calculate prediction metrics#
# input dataframe with two columns: predicted, original
Cell_type_stats_mean=function(DF_compare) {
  
  prediction_results <-table(DF_compare$predicted,DF_compare$original)
  # Define the true classes and predicted classes
  true_classes <- colnames(prediction_results)
  predicted_classes <- rownames(prediction_results)
  
  # Calculate true positives, false positives, false negatives, and true negatives
  TP <- numeric(length(true_classes))
  FP <- numeric(length(true_classes))
  FN <- numeric(length(true_classes))
  TN <- numeric(length(true_classes))
  
  for (i in 1:length(true_classes)) {
    true_class <- true_classes[i]
    
    true_class_index <- match(true_class, predicted_classes)
    
    TP[i] <- prediction_results[true_class_index, i]
    FN[i] <- sum(prediction_results[, i]) - TP[i]
    FP[i] <- sum(prediction_results[true_class_index, ]) - TP[i]
    TN[i] <- sum(prediction_results) - TP[i] - FP[i] - FN[i]
  }
  
  
  recall=100*TP/(TP+FN)
  recall[is.na(recall)] <- 0
  recall=mean(recall)
  
  precision=100*TP/(TP+FP)
  precision[is.na(precision)] <- 0
  precision=mean(precision)
  
  npv=100*TN/(TN+FN)
  npv[is.na(npv)] <- 0
  npv=mean(npv)
  
  spec=100*TN/(TN+FP)
  spec[is.na(spec)] <- 0
  spec=mean(spec)
  
  accuracy=100*(TP+TN)/(TP+TN+FP+FN)
  accuracy[is.na(accuracy)] <- 0
  accuracy=mean(accuracy)
  
  df=data.frame(Sensitivity=round(recall,digits = 2),
                PPV=round(precision,digits = 2),
                NPV=round(npv,digits = 2),
                Specificity=round(spec,digits = 2),
                Accuracy= round(accuracy,digits = 2)
  )
  df$f1_score=round(2*(df$PPV*df$Sensitivity)/(df$PPV+df$Sensitivity),digits=2)
  df
  
}

# Create a dataframe from the predicted cluster labels and original cluster labels
# Call Cell_type_stats_mean to calculate prediction metrics

harmonize_pred_stat_mean=function(new_clusters,orig){
  df_filter= data.frame(original=orig,predicted=new_clusters)
  df_filter$original=unfactor(df_filter$original)
  df_filter$predicted=unfactor(df_filter$predicted)
  print(sort(unique(df_filter$original)))
  print(sort(unique(df_filter$predicted)))
  
  data.frame(Cell_type_stats_mean(df_filter))
}

##### Cell type harmonization ##### 
# Cell type harmonization for sctypecelltypist library
sctype_harmonization=function(sctype_cells){
case_when(
    grepl("Monocytes",sctype_cells,ignore.case = T) ~"Mono",
    grepl("CD4+ T cells",sctype_cells,fixed = T) ~"T CD4",
    grepl("CD8+ T cells",sctype_cells,fixed = T) ~"T CD8",
    grepl("Dendritic cells",sctype_cells,fixed = T) ~"DC",
    grepl("B cells",sctype_cells,fixed = T) ~"B",
    grepl("Natural killer  cells",sctype_cells,fixed = T) ~"NK",
    grepl("HSC/MPP cells",sctype_cells,fixed = T) ~"HSC",
    grepl("CD4+ NKT-like cells",sctype_cells,fixed = T) ~"T unconv",
    grepl("CD8+ NKT-like cells",sctype_cells,fixed = T) ~"T unconv",
    grepl("γδ-T cells",sctype_cells,fixed = T) ~"T unconv",
    
    T~sctype_cells
  )
}
  
  