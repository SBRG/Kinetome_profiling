write.pdf <- FALSE

library(scales)

# load result created in "analyses/kapp_analysis/MOMENT_validation/MOMENT_ME_abu_preds_intersect_genes_Schmidt.R"
load("analyses/kapp_analysis/MOMENT_validation/MOMENT_saved_preds/Moment_Schmidt_data_exclude_transp_TRUE_intersect_4x_cutoff.rdata")

# load barplot_mechmod_perf
source("analyses/auxilliary_scripts/barplot_mech_model_perf.R")
# combine ggplots for MOMENT and ME validation

# also determines the order in plot
mods_to_to_plot <- c( "kcat_iv_ML_per_AS_per_s_NULL_med"      ,      "kcat_iv_ML_per_AS_per_s_repl"      ,
                      "kappmax_davidi_per_pp_per_s_NULL_med"   ,     "kappmax_davidi_per_pp_per_s_repl"  ,    
                      "kappmax_KO_ALE_per_pp_per_s_NULL_med"    ,    "kappmax_KO_ALE_per_pp_per_s_repl" ,   
                      "kappmax_KO_ALE_davidi_per_pp_per_s_NULL_med" ,"kappmax_KO_ALE_davidi_per_pp_per_s_repl"  
                      )

rename_map <- c( "kappmax_davidi_per_pp_per_s_repl"   = "Davidi et al. ensemble model" ,
                 "kappmax_davidi_per_pp_per_s_NULL_med" = "Davidi et al. median imputed"   ,  
                 "kappmax_KO_ALE_davidi_per_pp_per_s_repl" = "Davidi et al. + ALE KO ensemble model",   
                 "kappmax_KO_ALE_davidi_per_pp_per_s_NULL_med" = "Davidi et al. + ALE KO median imputed", 
                 "kappmax_KO_ALE_per_pp_per_s_repl" = "ALE KO ensemble model",
                 "kappmax_KO_ALE_per_pp_per_s_NULL_med" = "ALE KO median imputed",       
                 "kcat_iv_ML_per_AS_per_s_repl"      ="kcat in vitro ensemble model"   ,     
                 "kcat_iv_ML_per_AS_per_s_NULL_med"    = "kcat in vitro medain imputed"
)

head(intersected_perf$out_df)


fill_cols <- c( brewer_pal(palette = "Paired")(length(unique(intersected_perf$out_df$source)) ))

desat <- function(cols, sat=0.5) {
  X <- diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
  hsv(X[1,], X[2,], X[3,])
}

fill_cols <- desat(fill_cols , 0.65 )

legend_title <- ""#expression(bold("Source of "*italic("k")["eff"]*":")) #
legend_labels <- expression( atop(italic("k")[cat]* italic(" in vitro"),"median imputed") ,
                    atop(italic("k")[cat]* italic(" in vitro"),"ensemble model") ,
                    atop(italic("k")["app,max"],"Davidi et al. median imputed"),
                    atop(italic("k")["app,max"],"Davidi et al. ensemble model"),
                    atop(italic("k")["app,max"],"KO ALE median imputed"),
                    atop(italic("k")["app,max"],"KO ALE ensemble model"),
                    atop(italic("k")["app,max"],"Davidi et al. and KO ALE median imputed"),
                    atop(italic("k")["app,max"],"Davidi et al. and KO ALE ensemble model")
                    
              )

bp <- barplot_mechmod_perf(performance_df = intersected_perf$out_df ,
                     mods_to_to_plot = mods_to_to_plot,
                     rename_map = rename_map,
                     fill_cols = fill_cols ,
                     legend_labels = legend_labels , 
                     legend_title = legend_title
                     )

if (write.pdf){
  ggsave(bp$arranged_grob , width = 18, height = 18,units = c( "cm" ),
       device = "pdf" , 
       file = "Manuscript/Figures/MOMENT_ME_barplot/MOMENT_ME_barplot.pdf")
}

############################################################
# performance summaries for main text:

# kappmax vs in vitro

KO_kappmax <- intersected_perf$out_df %>%  filter(source == "kappmax_KO_ALE_per_pp_per_s_repl") %>% 
  #group_by( mech_model ) %>% 
  summarise( m.R2 = mean(R2),
             med.R2 = median(R2),
             m.RMSE = mean(RMSE),
             med.RMSE = median(RMSE)
             )


iv <- intersected_perf$out_df %>%  filter(source == "kcat_iv_ML_per_AS_per_s_repl") %>% 
  #group_by( mech_model ) %>% 
  summarise( m.R2 = mean(R2),
             med.R2 = median(R2),
             m.RMSE = mean(RMSE),
             med.RMSE = median(RMSE)
  )

KO_kappmax / iv

( iv - KO_kappmax ) / iv

##############################
# ML vs non-ML

kapp_max_null <- intersected_perf$out_df %>%  filter( grepl("kappmax.+NULL" , source ) ) %>% 
  #group_by( mech_model ) %>% 
  summarise( m.R2 = mean(R2),
             med.R2 = median(R2),
             m.RMSE = mean(RMSE),
             med.RMSE = median(RMSE)
  )

kapp_max_ML <- intersected_perf$out_df %>%  filter( grepl("kappmax.+_repl" , source ) ) %>% 
  #group_by( mech_model ) %>% 
  summarise( m.R2 = mean(R2),
             med.R2 = median(R2),
             m.RMSE = mean(RMSE),
             med.RMSE = median(RMSE)
  )


(kapp_max_null - kapp_max_ML) / kapp_max_null

# in vitro

iv_null <- intersected_perf$out_df %>%  filter( grepl("kcat_iv.+NULL" , source ) ) %>% 
  #group_by( mech_model ) %>% 
  summarise( m.R2 = mean(R2),
             med.R2 = median(R2),
             m.RMSE = mean(RMSE),
             med.RMSE = median(RMSE)
  )

iv_ML <- intersected_perf$out_df %>%  filter( grepl("kcat_iv.+_repl" , source ) ) %>% 
  #group_by( mech_model ) %>% 
  summarise( m.R2 = mean(R2),
             med.R2 = median(R2),
             m.RMSE = mean(RMSE),
             med.RMSE = median(RMSE)
  )


(iv_null - iv_ML) / iv_null


##############################
# KO vs Davidi

KO_kappmax <- intersected_perf$out_df %>%  filter(source == "kappmax_KO_ALE_per_pp_per_s_repl") %>% 
  group_by( mech_model ) %>% 
  summarise( m.R2 = mean(R2),
             med.R2 = median(R2),
             m.RMSE = mean(RMSE),
             med.RMSE = median(RMSE)
  ) 

davidi_kappmax <- intersected_perf$out_df %>%  filter(source == "kappmax_davidi_per_pp_per_s_repl") %>% 
  group_by( mech_model ) %>% 
  summarise( m.R2 = mean(R2),
             med.R2 = median(R2),
             m.RMSE = mean(RMSE),
             med.RMSE = median(RMSE)
  )

KO_kappmax
davidi_kappmax


(KO_kappmax[,2:5] - davidi_kappmax[,2:5]) / davidi_kappmax[,2:5]

