load("analyses/kapp_analysis/MOMENT_validation/MOMENT_saved_preds/Moment_Schmidt_data_exclude_transp_TRUE_intersect_4x_cutoff.rdata")

# load scatterplot_mechmod_perf
source("analyses/auxilliary_scripts/scatterplot_mech_model_perf.R")

conditions_to_plot <- c("Acetate", "Glucose")

performance_df <- intersected_perf$out_df 
performance_df <- performance_df[ !grepl("NULL",performance_df$source) ,  ] 

# adjust positions
performance_df$pos_x <- -1.8 
performance_df$pos_y <- -6.25

rename_map <- c( "kappmax_KO_ALE_per_pp_per_s_repl" = "kapp,max\nALE KO\nensemble model",
                 "kappmax_KO_ALE_per_pp_per_s_NULL_med" = "ALE KO\nmedian imputed",   
                 "kappmax_davidi_per_pp_per_s_repl"   = "kapp,max\nDavidi et al.\nensemble model" ,
                 "kappmax_davidi_per_pp_per_s_NULL_med" = "Davidi et al.\nmedian imputed"   ,  
                 "kappmax_KO_ALE_davidi_per_pp_per_s_repl" = "kapp,max\nDavidi et al. + ALE KO\nensemble model",   
                 "kappmax_KO_ALE_davidi_per_pp_per_s_NULL_med" = "Davidi et al. + ALE KO\nmedian imputed", 
                 "kcat_iv_ML_per_AS_per_s_repl"      ="kcat in vitro\nensemble model"   ,     
                 "kcat_iv_ML_per_AS_per_s_NULL_med"    = "kcat in vitro\nmedain imputed"
)


plot_out <- scatterplot_mechmod_perf( performance_df = performance_df, 
                          ab_df_overlapped = intersected_perf$ab_df_overlapped,
                          rename_map = rename_map,
                          conditions_to_plot = conditions_to_plot
                          )

plot_out

ggsave(plot_out, 
       file = "Manuscript/Figures/SI/MOMENT_ME_scatterplots/MOMENT_ME_scatterplots.png",
       width = 10, height =10
       )

