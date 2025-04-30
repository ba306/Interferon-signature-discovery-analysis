library(varhandle)
library(foreach)
library(ggpubr)
library(dplyr)
library("ggsci")
library(gridExtra)
library(cowplot)
library(reshape2)
library(rstatix)


text_size=40
signif_size=12
bracket_size=1.5

# Boxplots of mean signature scores of our IFN signaturs
# In discovery and validation datasets
# add t test stats
# for datasets with 2 samples each e.g. Rai and Fujiwara no stats added

data_dir="data/IFN_bulk/" 
save_dir="discovery/" 

# Load discovery dataset expression matrices + metadata
load(
  file = paste0(data_dir,
                "discovery_exp_meta_rdata.RData"))

# Load our IFN signatures
IFN_list_ind=readRDS(paste0(save_dir,"IFN_Aybey_finallist.RDS"))

# Metadata formatting for plotting 
fuji_df$rows=rownames(fuji_df)
colli_df$rows=rownames(colli_df)
colli_df$time=as.numeric(colli_df$time)
rai_df$rows=rownames(rai_df)
rai_df$time=as.numeric(rai_df$time)

# Function to calculate zscaled mean expression scores for a given signature
# given dataset + metadata

mean_meta_melt=function(dataset_name,test_data,
                        IFN_list,meta){
  
  data_j= t(scale(t(test_data), center=TRUE, scale=TRUE))
  score_list= IFN_list
  score_list=lapply(score_list, function(x){
    intersect(rownames(data_j),x)
  })
  test_dt=foreach(i=seq_along(score_list),.combine=cbind) %do% {
    ls=score_list[[i]]
    data.frame(post=
                 colMeans(data_j[ls,])
    )
  }
  colnames(test_dt)=names(score_list)
  
  if(is.factor(meta$Stim)==T){
    meta$Stim=unfactor(meta$Stim)}else{meta$Stim=meta$Stim}
  test_dt$rows=rownames(test_dt)
  
  
  test_dt_m=melt(test_dt)
  test_dt_m=left_join(test_dt_m,meta)
  test_dt_list=list(test_dt_m)
  names(test_dt_list)=dataset_name
  test_dt_list
}

#### Calculate scores for discovery datasets ####

Discovery_list=c(
  mean_meta_melt(dataset_name = "Jankowski",test_data = jank_c,
                 IFN_list = IFN_list_ind,meta = jank_df),
  mean_meta_melt(dataset_name = "Ziegler",test_data = ziegler_c,
                 IFN_list = IFN_list_ind,meta = ziegler_df),
  mean_meta_melt(dataset_name = "Fuji",test_data = fuji_c,
                 IFN_list = IFN_list_ind,meta = fuji_df),
  mean_meta_melt(dataset_name = "Colli",test_data = colli_c,
                 IFN_list = IFN_list_ind,meta = colli_df),
  mean_meta_melt(dataset_name = "Rai",test_data = rai_c,
                 IFN_list = IFN_list_ind,meta = rai_df))

Discovery_list=lapply(Discovery_list,function(d){
  d$Stim=gsub("untreated","Control",x = d$Stim)
  colnames(d)[2]="score"
  d
})


# function to calculate p values
add_stats=function(df){
  df_stats <- df %>%
    group_by(Stim) %>%
    t_test(value ~ score,p.adjust.method = "none")
  df_stats %>%
    add_xy_position(x = "Stim", dodge = 0.8)
  
}

#### Jankowski ####

plot_Jankowski=ggplot(data = Discovery_list$Jankowski,
                      aes(Stim,value,color=score))+
  geom_boxplot(outlier.colour=NULL)+
  # geom_point(position=position_dodge(0.5))+
  scale_color_nejm()+  theme_bw()+ theme(legend.position="top",
                                         text = element_text(size = text_size))+
  labs(title = "Jankowski",color="IFN signatures")+
  xlab("Stimulation")+rremove("ylab")+ 
  stat_pvalue_manual(
    add_stats(Discovery_list$Jankowski), label = "p.adj.signif", tip.length = 0.01,hide.ns = T,size=signif_size,
    bracket.size = bracket_size
  )

#### Ziegler donor 1 ######
dose.labs <- paste0(unique(Discovery_list$Ziegler$Dose)," ng/ml")
names(dose.labs) <-   unique(Discovery_list$Ziegler$Dose)

ziegler_df=Discovery_list$Ziegler
add_untreated_a=ziegler_df[ziegler_df$Stim=="Control",]
add_untreated_a$Stim=gsub("Control","IFNa",add_untreated_a$Stim)

add_untreated_g=ziegler_df[ziegler_df$Stim=="Control",]
add_untreated_g$Stim=gsub("Control","IFNg",add_untreated_g$Stim)

ziegler_df2=ziegler_df[ziegler_df$Stim!="Control",]
ziegler_df2=rbind(ziegler_df2,add_untreated_a,add_untreated_g)

# add stats

zieglerdonor1_stats <- ziegler_df2 %>%
  group_by(Dose,Stim) %>%
  t_test(value ~ score,p.adjust.method = "none")

zieglerdonor1_stats=zieglerdonor1_stats %>%
  add_xy_position(x = "Dose", dodge = 0.8)

plot_Ziegler_donor1=ggboxplot(ziegler_df2,x="Dose",y="value",
                              width = 0.4  ,color="score",
                              outlier.shape = NULL)+
  scale_color_nejm()+ theme_bw()+ theme(legend.position = "top", text = element_text(size = text_size),
                                        strip.text.y = element_text(size = text_size),
                                        strip.text.x = element_text(size = text_size)) + 
  labs(title = "Ziegler donor 1",color="IFN signatures")+
  xlab("Dose [ng/ml]")+
  # ylab("Mean signature score")+
  facet_wrap(.~Stim,nrow = 1)+ rremove("ylab")+
  stat_pvalue_manual(
    zieglerdonor1_stats, label = "p.adj.signif", tip.length = 0.01,hide.ns = T,size=signif_size,
    bracket.size = bracket_size
  )


#### Fujiwara #####

# didnt add stats since two samples each

plot_Fuji=ggplot(data = Discovery_list$Fuji,
                 aes(Stim,value,color=score))+
  geom_boxplot(outlier.colour=NULL)+
  # geom_point(position=position_dodge(0.7),size=1)+
  scale_color_nejm()+  theme_bw()+  theme(legend.position="top",text = element_text(size = text_size))+
  labs(title = "Fujiwara",color="IFN signatures")+
  xlab("Stimulation")+rremove("ylab") 
# stat_pvalue_manual(
#   add_stats(Discovery_list$Fuji), label = "p.adj.signif", tip.length = 0.01,hide.ns = T
# )



#### Colli ####
time.labs <- paste0(unique(Discovery_list$Colli$time)," h")
names(time.labs) <-   unique(Discovery_list$Colli$time)


# add stats

colli_stats <-  Discovery_list$Colli %>%
  group_by(time,Stim) %>%
  t_test(value ~ score,p.adjust.method = "none")

colli_stats=colli_stats %>%
  add_xy_position(x = "Stim", dodge = 0.8)

plot_colli=ggplot(data = Discovery_list$Colli,
                  aes(Stim,value,color=score))+
  geom_boxplot(outlier.colour=NULL)+
  # geom_point(position=position_dodge(0.7),size=1)+
  scale_color_nejm()+  theme_bw()+ theme(legend.position="top",text = element_text(size = text_size),
                                         strip.text.y = element_text(size = text_size),
                                         strip.text.x = element_text(size = text_size))+
  labs(title = "Colli",color="IFN signatures")+
  xlab("Stimulation")+
  # ylab("Mean signature score")+
  facet_wrap(.~time,scales = "free_x",
             labeller = labeller(time = time.labs))+ rremove("ylab")+ 
  stat_pvalue_manual(
    colli_stats, label = "p.adj.signif", tip.length = 0.01,hide.ns = T,size=signif_size,
    bracket.size = bracket_size
  )


#### Rai #####
time.labs <- paste0(unique(Discovery_list$Rai$time)," h")
names(time.labs) <-   unique(Discovery_list$Rai$time)

# didnt add stats since two samples each time point 

plot_Rai=ggplot(data = Discovery_list$Rai,
                aes(Stim,value,color=score))+
  geom_boxplot(outlier.colour=NULL)+
  # geom_point(position=position_dodge(0.7),size=1)+
  scale_color_nejm()+  theme_bw()+  theme(legend.position="top",text = element_text(size = text_size))+
  labs(title = "Rai",color="IFN signatures")+
  xlab("Stimulation")+
  # ylab("Mean signature score")+
  facet_wrap(.~time,scales = "free_x",
             labeller = labeller(time = time.labs))+ rremove("ylab")
# Get legend 
# legend <- cowplot::get_legend(
#   plot_Fuji
# )

legend=get_plot_component(plot_Rai, 'guide-box', return_all = TRUE)
cowplot::ggdraw(legend)

# Upper plots for discovery

upper_discovery=
  plot_grid(
    plot_Jankowski+ theme(legend.position="none"),
    plot_Ziegler_donor1+ theme(legend.position="none"),
    # legend,
    # labels = c("A", "B","C"),
    nrow = 1,
    rel_widths = c(0.7, 2)
  )

# Bottom plots for discovery

bottom_discovery=
  plot_grid(
    plot_Fuji+ theme(legend.position="none"),
    plot_colli+ theme(legend.position="none"),
    plot_Rai+ theme(legend.position="none"),
    # labels = c("C", "D","E"), 
    nrow = 1,
    rel_widths = c(0.7,1,1)
  )


#### Validation datasets ####
# Load datasets expression + metadata

load(
  file = paste0(data_dir,
                "validation_exp_meta_rdata.RData"))

#### Calculate scores for validation datasets ####

Test_list=c(
  mean_meta_melt(dataset_name = "ziegler_donor_2",test_data = ziegler_donor2_list[[1]],
                 IFN_list = IFN_list_ind,meta = ziegler_donor2_list[[2]]),
  mean_meta_melt(dataset_name = "ziegler_beas",test_data = ziegler_beas_list[[1]],
                 IFN_list = IFN_list_ind,meta = ziegler_beas_list[[2]]),
  mean_meta_melt(dataset_name = "lee",test_data = lee_list[[1]],
                 IFN_list = IFN_list_ind,meta = lee_list[[2]]),
  mean_meta_melt(dataset_name = "devlin",test_data = devlin_list[[1]],
                 IFN_list = IFN_list_ind,meta = devlin_list[[2]]))

Test_list=lapply(Test_list,function(d){
  d$Stim=gsub("untreated","Control",x = d$Stim,ignore.case = T)
  d$Stim=gsub("Buffer","Control",x = d$Stim,ignore.case = T)
  d$Stim=gsub("IFNA","IFNa",x = d$Stim)
  d$Stim=gsub("IFNG","IFNg",x = d$Stim)
  d$Stim=gsub("IFNgamma","IFNg",x = d$Stim)
  d$Stim=gsub("IFNbeta","IFNb",x = d$Stim)
  d$Stim=gsub("IFNl3","IFNlambda",x = d$Stim)
  colnames(d)[2]="score"
  d
})

#### Ziegler donor 2 ####

ziegler_df_donor2=Test_list$ziegler_donor_2
add_untreated_a=ziegler_df_donor2[ziegler_df_donor2$Stim=="Control",]
add_untreated_a$Stim=gsub("Control","IFNa",add_untreated_a$Stim)

add_untreated_g=ziegler_df_donor2[ziegler_df_donor2$Stim=="Control",]
add_untreated_g$Stim=gsub("Control","IFNg",add_untreated_g$Stim)

ziegler_df2_donor2=ziegler_df_donor2[ziegler_df_donor2$Stim!="Control",]
ziegler_df2_donor2=rbind(ziegler_df2_donor2,add_untreated_a,add_untreated_g)


# add stats
zieglerdonor2_stats <- ziegler_df2_donor2 %>%
  group_by(Dose,Stim) %>%
  t_test(value ~ score,p.adjust.method = "none")

zieglerdonor2_stats=zieglerdonor2_stats %>%
  add_xy_position(x = "Dose", dodge = 0.8)

plot_Ziegler_donor2=ggboxplot(ziegler_df2_donor2,x="Dose",y="value",
                              width = 0.4  ,color="score",
                              facet.by = "Stim",outlier.shape = NULL)+
  scale_color_nejm()+  theme_bw()+
  labs(title = "Ziegler donor 2",color="IFN signatures")+ theme(legend.position="top",
                                                                text = element_text(size = text_size),
                                                                strip.text.y = element_text(size = text_size),
                                                                strip.text.x = element_text(size = text_size))+
  xlab("Dose [ng/ml]")+rremove("ylab")+ 
  stat_pvalue_manual(
    zieglerdonor2_stats, label = "p.adj.signif", tip.length = 0.01,hide.ns = T,size=signif_size,
    bracket.size = bracket_size
  )

#### Ziegler Beas cell line ###

ziegler_df_beas=Test_list$ziegler_beas
add_untreated_a=ziegler_df_beas[ziegler_df_beas$Stim=="Control",]
add_untreated_a$Stim=gsub("Control","IFNa",add_untreated_a$Stim)

add_untreated_g=ziegler_df_beas[ziegler_df_beas$Stim=="Control",]
add_untreated_g$Stim=gsub("Control","IFNg",add_untreated_g$Stim)

ziegler_df2_beas=ziegler_df_beas[ziegler_df_beas$Stim!="Control",]
ziegler_df2_beas=rbind(ziegler_df2_beas,add_untreated_a,add_untreated_g)


# add stats
zieglerbeas_stats <- ziegler_df2_beas %>%
  group_by(Dose,Stim) %>%
  t_test(value ~ score,p.adjust.method = "none")

zieglerbeas_stats=zieglerbeas_stats %>%
  add_xy_position(x = "Dose", dodge = 0.8)

plot_Ziegler_beas=ggboxplot(ziegler_df2_beas,x="Dose",y="value",
                            width = 0.4  ,color="score",
                            facet.by = "Stim",outlier.shape = NULL)+
  scale_color_nejm()+  theme_bw()+theme(legend.position="top",text = element_text(size = text_size),
                                        strip.text.y = element_text(size = text_size),
                                        strip.text.x = element_text(size = text_size))+
  labs(title = "Ziegler BEAS-2B",color="IFN signatures")+
  xlab("Dose [ng/ml]")+rremove("ylab")+ 
  stat_pvalue_manual(
    zieglerbeas_stats, label = "p.adj.signif", tip.length = 0.01,hide.ns = T,size=signif_size,
    bracket.size = bracket_size
  )


### Lee ####

plot_lee=ggplot(data = Test_list$lee,
                aes(Stim,value,color=score))+
  geom_boxplot(outlier.colour=NULL)+
  # geom_point(position=position_dodge(0.7),size=1)+
  scale_color_nejm()+  theme_bw()+theme(legend.position="top",text = element_text(size = text_size))+
  labs(title = "Lee",color="IFN signatures")+
  xlab("Stimulation")+rremove("ylab")+ 
  stat_pvalue_manual(
    add_stats(Test_list$lee), label = "p.adj.signif", tip.length = 0.01,hide.ns = T,size=signif_size,
    bracket.size = bracket_size
  )

#### Devlin ####

plot_devlin=ggplot(data = Test_list$devlin,
                   aes(Stim,value,color=score))+
  geom_boxplot(outlier.colour=NULL)+
  # geom_point(position=position_dodge(0.7),size=1)+
  scale_color_nejm()+  theme_bw()+theme(legend.position="top",text = element_text(size = text_size))+
  labs(title = "Devlin",color="IFN signatures")+
  xlab("Stimulation")+
  rremove("ylab")+ 
  stat_pvalue_manual(
    add_stats(Test_list$devlin), label = "p.adj.signif", tip.length = 0.01,hide.ns = T,size=signif_size,
    bracket.size = bracket_size
  )

### All validation plots ####

test_plots=plot_grid(
  plot_Ziegler_donor2+ theme(legend.position="none"),
  plot_lee+ theme(legend.position="none"),
  plot_Ziegler_beas+ theme(legend.position="none"),
  plot_devlin+ theme(legend.position="none"),
  # labels = "AUTO", nrow = 2,
  rel_widths = c(2,1,2,1)
)

#### Put discovery and validation boxplots together ####

pdf(paste0(save_dir,"Discovery_validation_boxplots_stats.pdf"),
    width = 30,
    height = 35)

all_discovery=plot_grid(upper_discovery, 
                        bottom_discovery,ncol = 1,
                        rel_heights = c(2,2))

# Add a common y-axis label
common_label <- ggdraw() +
  draw_label("Mean signature score", size = 40, angle = 90,fontface = "plain")

# Combine the grid and common label
all_discovery_common=plot_grid(common_label,all_discovery,ncol = 2, align = "v",rel_widths = c(0.02,1))
test_plots_common=plot_grid(common_label,test_plots,ncol = 2, align = "v",rel_widths = c(0.02,1))

print(plot_grid(legend,
                all_discovery_common,
                ggplot() + theme_void(),  # Empty ggplot for space without gray area
                test_plots_common,
                ncol = 1,
                rel_heights = c(0.1, 2, 0.2, 2),  # Adjust spacer height as needed
                labels = c("", "A", "", "B"),
                scale = 1,
                vjust = 0.1,
                label_size = 60))


dev.off()



