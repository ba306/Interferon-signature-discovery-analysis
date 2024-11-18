#coherence score
coherence_score_val <- function(genes, pearson_cor) {
  cor <- 0
  missing_genes <- vector()
  for (gene1 in 1:length(genes)) {
    for (gene2 in 1:length(genes)) {
      if (gene1 < gene2) {
        if (!(genes[gene1] %in% row.names(pearson_cor))) {
          missing_genes <- c(missing_genes, genes[gene1])
          next
        } else if (!(genes[gene2] %in% row.names(pearson_cor))) {
          missing_genes <- c(missing_genes, genes[gene2])
          next
        } else {
          if (!is.na(pearson_cor[genes[gene1], genes[gene2]])) {
            cor <- cor + pearson_cor[genes[gene1], genes[gene2]] 
          }
        }
      }
    }
  }
  gene_length <- length(genes) - length(unique(missing_genes) )
  cs_score <- 2*cor / (gene_length * (gene_length-1)) 
  if (is.na(cs_score)) {
    cs_score <- 0
  }
  # print(genes)
  # print(cs_score)
  # print(unique(missing_genes))
  return (cs_score)
}
