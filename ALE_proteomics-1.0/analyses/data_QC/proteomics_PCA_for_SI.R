########
# WARNING run on fresh session, plyr/ dplyr conflicts can occur

####
# parameters
# use top 75 percent most abundant prots
abundance_cutoff <- 0.25 # percentile/100 for under which smaples are discarded
write_plots <- FALSE

library(plyr);library(dplyr)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(scales)
library(Rtsne)
library(gridExtra)

# created in "compute_abundance_proxies.R"
load("data/proteomics/ALE_Proteomics_1_26_2019_ab_proxies.rdata")
source("data/proteomics/sample_id_name_map.R")

##########################################################
# sample correlation in top3

# use top 75 percent most abundant prots ()
top3_quant_cutoff <- quantile( ab_proxies$top3 , abundance_cutoff )
log10(top3_quant_cutoff)

# remove UPS for this part of the analysis
top3 <- ab_proxies %>% #select( - ibaq ) %>% 
  filter( !grepl("ups$", Master.Protein.Accessions ) ) %>% 
  filter(sample_id %in% names(sample_id_map)) # focus on samples used in the kappmax study

# filter low-abundance proteins
top3_f <- top3 %>% filter( top3 > top3_quant_cutoff )  %>% filter( !grepl("ups$", Master.Protein.Accessions ) )
top3_f_repmean <- top3_f %>% group_by(sample_id , Master.Protein.Accessions) %>% summarise(top3_repmean = mean(top3,na.rm = TRUE)) 

top3_wide <- acast( top3 , formula =  Master.Protein.Accessions ~ sample_id + bio_rep , value.var = "top3")
top3_f_wide <- acast( top3_f , formula =  Master.Protein.Accessions ~ sample_id + bio_rep , value.var = "top3")
top3_f_repmean_wide <- acast( top3_f_repmean , formula =  Master.Protein.Accessions ~ sample_id , value.var = "top3_repmean")
#top3_wide <- acast( top3 , formula =  Master.Protein.Accessions ~ sample_id , value.var = "top3")

# only use complete obs
top3_f_wide_compl_obs <- top3_f_wide[ apply(top3_f_wide,1,function(x){!any(is.na(x))}), ]
top3_f_repmean_wide_compl_obs <- top3_f_repmean_wide[ apply(top3_f_repmean_wide,1,function(x){!any(is.na(x))}), ]

########################################################
# PCA
library(ggbiplot)
library(ggthemes)

top3_f_wide_pca <- top3_f_repmean_wide_compl_obs
dim(top3_f_wide_pca)

colnames(top3_f_wide_pca) <- sample_id_map[ match( colnames(top3_f_wide_pca), names(sample_id_map) ) ]

condition_groups <- sub( " \\d$" ,  "" , colnames(top3_f_wide_pca))

pc <- prcomp( t( scale(log10( top3_f_wide_pca + 1 )) ) , center = TRUE , scale. = TRUE )

pca_biplot <- ggbiplot(pc, choices = 1:2 , groups = condition_groups , labels =  colnames(top3_f_wide_pca) , labels.size = 3, alpha = 0.5 , 
                       ellipse = FALSE, circle = FALSE, varname.size=0.1, var.axes = FALSE)
pca_biplot <- pca_biplot + ggthemes::theme_base() +
  guides(color = guide_legend(title = "background of\nevolved strain") ) +
  theme(plot.background=element_rect( colour=NA))

pca_biplot

if(write_plots){
  pdf("Manuscript/Figures/SI/proteomics_PCA/proteomics_PCA.pdf")
  pca_biplot
  dev.off()
  
  ggsave(pca_biplot , filename = "Manuscript/Figures/SI/proteomics_PCA/proteomics_PCA.png" )
}

# most important genes in PC1
pc1_ord <- pc$rotation[,"PC1"][order(abs(pc$rotation[,"PC1"]),decreasing = T )]
hist(pc1_ord)
pc4_ord <- pc$rotation[,"PC4"][order(abs(pc$rotation[,"PC4"]),decreasing = T )]

#########################################################
set.seed(162)  
tsne <- Rtsne( t( scale(log10(top3_f_wide_pca+1) )) , perplexity = 5  ) #, check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.5, dims=2)

ggplot( as.data.frame(tsne$Y), aes(x=V1, y=V2 ,col = condition_groups )) +  
  geom_point() + geom_text_repel( aes( label = colnames(top3_f_wide_pca)  ) , size = 4 ) + 
  xlab("t-SNE dimension 1") + ylab("t-SNE dimension 2") + 
  ggthemes::theme_base() + theme(legend.position="none")


#########################################################
# PCAs: BOP + individual conditions
ind_biplot <- vector(length = length(unique(condition_groups)) , mode = "list")
names(ind_biplot) <- unique(condition_groups)
for ( i in seq_along(unique(condition_groups) ) ){
  col_log <-  condition_groups %in% c("BOP",unique(condition_groups)[i])  
  ind_biplot[[i]] <- ggbiplot(prcomp( t( scale(log10( top3_f_wide_pca[ ,col_log] + 1 ) )) , center = TRUE , scale. = TRUE ), 
                              choices = 1:2 , groups = condition_groups[ col_log  ] , labels =  colnames(top3_f_wide_pca)[col_log] , labels.size = 3, alpha = 0.5 , 
                              ellipse = FALSE, circle = FALSE, varname.size=0.1, var.axes = FALSE) + ggthemes::theme_base()
  
}

do.call("grid.arrange", c(ind_biplot[names(ind_biplot)!="BOP"] , ncol =4) )
#########################################################
# densities

top3_f %>% mutate( sample_id_short = substr( sample_id ,1,3) , ALE_num = substr( sample_id ,nchar(sample_id),nchar(sample_id)) ) %>%  #filter(sample_id == "BOP122") %>%  
  ggplot(aes(x = log10(top3) , col = ALE_num , linetype = bio_rep) ) + 
  geom_density() + 
  #stat_bin(geom="step", position="identity") +
  facet_wrap( . ~ sample_id_short)

# note outliers:
top3_f %>% ggplot(aes(x = log10(top3) , col = sample_id ) ) + geom_density() #+ geom_text( aes(label = "sample_id" , y = 0.5 ))

