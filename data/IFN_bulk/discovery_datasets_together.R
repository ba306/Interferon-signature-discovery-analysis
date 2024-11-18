# Combine meta and expression data from discovery datasets 
# log2 transform expression data 

save_dir="./discovery/"

data_dir="./data/IFN_bulk/" 

#### Read  Datasets  ####

# Read expression data
ziegler=readRDS(paste0(data_dir,"Ziegler/",
                       "ziegler_donor1.RDS"))[[1]] #raw

jank=readRDS(paste0(data_dir,"Jankowski/",
                    "jankowski.RDS"))[[1]] #raw

fuji=readRDS(paste0(data_dir,"Fujiwara/",
                    "fujiwara.RDS"))[[1]] # fpkm

colli=readRDS(paste0(data_dir,"Colli/",
                     "colli.RDS"))[[1]] # tpm

rai=readRDS(paste0(data_dir,"Rai/",
                   "rai.RDS"))[[1]] # raw


# Read metadata

ziegler_df=readRDS(paste0(data_dir,"Ziegler/",
                          "ziegler_donor1.RDS"))[[2]]

jank_df=readRDS(paste0(data_dir,"Jankowski/",
                       "jankowski.RDS"))[[2]] 

fuji_df=readRDS(paste0(data_dir,"Fujiwara/",
                       "fujiwara.RDS"))[[2]]

colli_df=readRDS(paste0(data_dir,"Colli/",
                        "colli.RDS"))[[2]]

rai_df=readRDS(paste0(data_dir,"Rai/",
                      "rai.RDS"))[[2]] 

# log2(x+1) transformation

fuji_c=log2(fuji+1)

ziegler_c=log2(ziegler+1)

jank_c=log2(jank+1)

colli_c=log2(colli+1)

rai_c=log2(rai+1)

#save for plotting purposes- boxplots
save(ziegler_c,ziegler_df,
     jank_c,jank_df,
     fuji_c,fuji_df,
     colli_c,colli_df,
     rai_c,rai_df,
     file = paste0(data_dir,
                   "discovery_exp_meta_rdata.RData"))
