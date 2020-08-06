write_pdf <- TRUE

# boxplot of CV model performances for in vivo and in vitro data


library(tidyverse)
library(ggthemes)
library(gridExtra)
library(egg)

# load "resamps", created in "get_ensemble_resamples.r"
load( "analyses/kapp_analysis/ML/saved_resamples/kappmax_preds_and_models_5000_dl_mods_3600_s_max_runtime_4x_cutoff_resamps_org_epochs.rdata")

resamps$Rsquared[ resamps$Rsquared < 0  ] <- 0

# remove the lm as it was not used in the ensemble:
resamps <- resamps[ ! grepl("_lm" , resamps$ML_model) , ]

table(resamps$ML_model)

keff_source_map <- expression( "kcat_iv_ML_per_AS_per_s" = atop(italic(k)[cat] ,"in vitro") ,
                      "kappmax_davidi_per_pp_per_s" =  atop(italic(k)["app,max"],"Davidi et al.")  , 
                      "kappmax_KO_ALE_per_pp_per_s" =  atop(italic(k)["app,max"],"KO ALE") ,
                      "kappmax_KO_ALE_davidi_per_pp_per_s" = atop(italic(k)["app,max"],atop("Davidi et al.","and KO ALE") )
                       )
keff_source_map_for_labels <- keff_source_map
names(keff_source_map_for_labels) <- NULL

resamps$keff_source <- keff_source_map[ match( resamps$keff_source , names(keff_source_map) ) ]


#set levels for bar order
resamps$keff_source <- factor(resamps$keff_source , levels = keff_source_map )


lapply( resamps , class )  

fill_cols <-  c( "#538DB4",  "#59A055", "#E36062", "#FFAC59" )


# aggregate across all models and show performance per input type

RMSE_agg <- resamps %>%
  ggplot( aes( x = keff_source , y = RMSE, fill = keff_source) ) +
  geom_boxplot(outlier.size=1,outlier.shape=4) +
  ylim(0,NA) +
  ggthemes::theme_few() + 
  scale_x_discrete( labels = rep(NULL,length(keff_source_map_for_labels))  ) + 
  scale_fill_manual(values = fill_cols) +
  theme(legend.position="none" ,  axis.title.x=element_blank() , 
        plot.margin = unit(c(0.1,0.2,0.1,0.1), "cm") ,
        axis.title.y = element_text( size = 10 ) 
        ) + 
  scale_y_continuous(labels=function(x) sprintf("%.2f", x)) + 
  ylab(expression("cross-validated RMSE (on log"[10]*" scale)"^" ")) 

R2_agg<- resamps %>%
  ggplot( aes( x = keff_source , y = Rsquared, fill = keff_source) ) +
  geom_boxplot(outlier.size=1,outlier.shape=4) +
  ylim(0,1) +
  ggthemes::theme_few() + 
  scale_x_discrete( labels = keff_source_map_for_labels  ) + 
  scale_fill_manual(values = fill_cols) +
  theme(legend.position="none",
        plot.margin = unit(c(0,0.2,1,0.1), "cm"),
        axis.title.x= element_text( size = 10 ,vjust = -2) ,
        axis.title.y = element_text( size = 10 ,vjust = 1) 
        ) +
  xlab( "source of turnover numbers" ) + 
  ylab(expression("cross-validated "*italic(R)^2*"(on log"[10]*" scale)") )

arranged_plot <- ggarrange(RMSE_agg,R2_agg  )


if(write_pdf){
  ggsave(file = "Manuscript/Figures/ML_CV_performance/ML_CV_performance.pdf" , 
         plot = arranged_plot , 
         width = 11, height = 17 , unit = "cm",
         device = "pdf")
}
