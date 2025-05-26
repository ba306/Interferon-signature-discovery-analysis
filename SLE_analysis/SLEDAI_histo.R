# Distribution of SLEDAI scores

pdf("./SLE_analysis/SLEDAI_histos.pdf",width = 5,height = 3)
par(mfrow = c(1, 3))  # 1 row, 3 columns
load("~/Interferon-signature-discovery-analysis/SLE_analysis/GSE121239_mtx_meta.RData")
meta_1=meta[grepl("v1", meta$title) & meta$`disease state:ch1`=="Systemic Lupus Erythematosus",]
meta_1$`sledai:ch1` %>% as.numeric() %>% hist(main="GSE121239") 

load("~/Interferon-signature-discovery-analysis/SLE_analysis/GSE49454_mtx_meta.RData")
meta_2=meta[grepl("V1", meta$title) & meta$`group:ch1`=="SLE",]

meta_2$`sledai:ch1` %>% as.numeric() %>% hist(main="GSE49454") 

load("~/Interferon-signature-discovery-analysis/SLE_analysis/GSE65391_mtx_meta.RData")
meta_3=meta[meta$`visit:ch1`==1 & meta$`disease state:ch1` =="SLE",]
meta_3$`sledai:ch1` %>% as.numeric() %>% hist(main="GSE65391") 

dev.off()


