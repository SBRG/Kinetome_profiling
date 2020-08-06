library(plyr); library(dplyr)
library(tidyverse)
library(reshape2)
library(mice)

write_plots <- FALSE

options(tibble.width = Inf,
        max.print = 1200
) # for printing tbls

# created in "calc_kapp_MFA.R"
t3_1tn_mfa <- read.csv("output_data/kappmax/KO_MFA_kappmax/KO_MFA_kappmax_full.csv" , stringsAsFactors = FALSE )

table(t3_1tn_mfa$sample_id)
source("data/proteomics/sample_id_name_map.R")

for(i in seq_along(sample_id_map) ){
  t3_1tn_mfa$sample_id[ names(sample_id_map)[i] == t3_1tn_mfa$sample_id ] <- sample_id_map[i]
  }

kapp_w <- t3_1tn_mfa %>% filter(!is.na(kapp_per_s)) %>% 
  group_by( sample_id , react_id ) %>%  # average replicates
  summarise( kapp_per_s = mean ( kapp_per_s, na.rm = TRUE) ) %>% 
  #mutate( sample_id_merge = paste0(sample_id," ",bio_rep)) %>% 
  select(sample_id , kapp_per_s , react_id)  %>% 
  acast( react_id ~  sample_id , value.var ="kapp_per_s")

kapp_w_comp <- kapp_w

# only use complete reactions
kapp_w_comp <- kapp_w[ apply(kapp_w,1, function(x){! any ( is.na(x) ) } ) , ]

tail(kapp_w)
tail(kapp_w_comp)

any(kapp_w_comp<=0)

#####################################
#  PCA

pca <-prcomp( t( log10(kapp_w_comp)) ,center = TRUE , scale = TRUE )

a <- predict(pca)

#library(ggbiplot)
source("analyses/auxilliary_scripts/ggbiplot_mod.r")

library(ggthemes)
library(ggrepel)

##############
# inspect PCs
pc1 <- pca$rotation[,1]
pc2 <- pca$rotation[,2]

sort(pc1)
sort(pc2)

hist(pc1)
hist(pc2)

##############
#plot 


condition_groups <- substr( colnames(kapp_w_comp) , 1 , nchar(colnames(kapp_w_comp) ) - 2)

ggbi <- ggbiplot(pca , 
         labels = colnames(kapp_w_comp), 
         groups = condition_groups,
         labels.size = 4.5, 
         alpha = 0.75 , 
         ellipse = FALSE, 
         circle = FALSE, 
         varname.size=0.1, 
         var.axes = FALSE #,
         #point_size = 4.5
         ) + 
  #ggplot2::theme_classic() +
  ggthemes::theme_base() +
  guides(color = guide_legend(title = "background of\nevolved strain") ) +
  theme(plot.background=element_rect( colour=NA))
  
ggbi

if(write_plots){
  pdf(file = "Manuscript/Figures/kapp_structure/kapp_biplot.pdf" , width = 7, height = 7 )
  ggbi
  dev.off()
}

##########################
# tsne
library(Rtsne)


set.seed(162)  
tsne <- Rtsne( t( scale(log10(kapp_w_comp) )) , perplexity = 5  ) 

ggplot( as.data.frame(tsne$Y), aes(x=V1, y=V2 ,col = condition_groups )) +  
  geom_point() + geom_text_repel( aes( label = colnames(kapp_w_comp) ) )

