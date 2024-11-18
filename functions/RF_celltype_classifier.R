# 1. Benchmarking functions/wrappers for RF training and prediction &
# 2. Prediction metrics calculation


library(randomForest)
library(varhandle)
library(foreach)
library(dplyr)
library(Seurat)
set.seed(42)

##### RF signature training functions ##### 

# Calculate and save z-scores of signature genes
zscores_func=function(data,assay,common_genes){
  data=data[common_genes,]
  data <- ScaleData(data, features = rownames(data))
  t(data@assays[[assay]]@scale.data)
  
}

# Training genes from the signatures using RF
RF_signatures=function(geneset, #name of the gene set, chr
                       signature_genes, #signature genes as in character vector or numeric (for random genes)
                       assay="RNA",
                       reference=NULL,
                       ref_col=NULL,
                       query=NULL,
                       ref_name=NULL,
                       query_name=NULL,
                       save_dir=NULL){
  if(is.numeric(signature_genes)){
    common_genes=intersect(
      rownames(query),
      rownames(reference)
    )
    
    signature_genes=sample(common_genes, signature_genes)
    
  }else{
    signature_genes=signature_genes}
  
  #  Test dataset
  print(geneset)
  genes_take_cell_types=intersect(rownames(query),signature_genes)
  print("Query dataset genes / Signature genes")
  print(paste(length(genes_take_cell_types),"/", length(unique(signature_genes))))
  query=query[genes_take_cell_types,]
  
  #Reference genes
  genes_take_cell_types=intersect(rownames(reference),signature_genes)
  print("Reference dataset genes / Signature genes")
  print(paste(length(genes_take_cell_types),"/", length(unique(signature_genes))))
  reference=reference[genes_take_cell_types,]
  
  # reference_seurat=reference_seurat[genes_take_cell_types,]
  common_genes=intersect(
    rownames(reference),
    rownames(query)
  )
  print(paste0("Common genes: ",length(common_genes)))
  
  
  hao_mean=zscores_func(reference,assay,common_genes)
  
  hao_mean=cbind(hao_mean,data.frame(cell_type=reference@meta.data[,ref_col]))
  
  hao_mean=data.frame(hao_mean)
  hao_mean$cell_type=factor(hao_mean$cell_type)
  
  # Seperate cells into train and test set 67:33
  ss <- sample(1:2,size=nrow(hao_mean),replace=TRUE,prob=c(0.67,0.33))
  
  train= hao_mean[ss==1,]
  test= hao_mean[ss==2,]
  
  # Train the data 
  start.time <- Sys.time()
  
  rf <- randomForest(
    cell_type ~ .,
    data=train,importance = TRUE
  )
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  
  
  file_extension=paste0("rf_",ref_name, "_", query_name, "_", geneset, "_", length(common_genes))
  
  saveRDS(rf,paste0(save_dir,file_extension,".RDS"))
  
  print("Training evaluation")
  print(harmonize_pred_stat_mean(rf$predicted,rf$y))
  
  # Predict cell types using trained RF model on test set
  pred = predict(rf, test)
  
  print("Training test data evaluation")
  print(harmonize_pred_stat_mean(pred,test$cell_type))
  
  #
  ko_mean=zscores_func(query,"RNA",common_genes)
  colnames(ko_mean)=gsub("-",".",colnames(ko_mean),fixed = T)
  #
  pred = predict(rf, ko_mean)
  pred=as.data.frame(pred)
  colnames(pred)="Predictions"
  
  
  write.csv(pred,file = 
              paste0(save_dir,file_extension,
                     ".csv")
  )
  
}


