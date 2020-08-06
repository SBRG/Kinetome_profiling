
scatterplot_mechmod_perf <- function(performance_df , # source , mech_model, condition, R2, RMSE
                                 #mods_to_to_plot, # which models with keff sources listed in "source"
                                 rename_map, # 
                                 ab_df_overlapped,
                                 #fill_cols, legend_labels , 
                                 #legend_title,
                                 conditions_to_plot  
) {
    library(ggplot2)
    library(ggthemes)
    library(gridExtra)
    library(scales)
    library(ggrepel)
    library(dplyr)

    ab_df_overlapped$source <- as.character(ab_df_overlapped$source)
    ab_df_overlapped$source <- factor( rename_map[ match( ab_df_overlapped$source , names(rename_map)) ] , levels = unique(rename_map))
  
    performance_df$source <- as.character(performance_df$source)
    performance_df$source <- factor( rename_map[ match( performance_df$source , names(rename_map)) ] , levels = unique(rename_map))
    
      
    ###########################################
  
    ab_overlapped <- merge( x = ab_df_overlapped ,  y = performance_df , 
                            by.x = c("source", "mech_model" , "growth_condition") ,
                            by.y = c("source", "mech_model" , "condition")  )
    


    performance_df_s <- performance_df %>% filter( condition %in% conditions_to_plot )
    names(performance_df_s)[names(performance_df_s) == "condition" ] <- "growth_condition"

    plot_out <- ab_overlapped %>% 
      filter( growth_condition %in% conditions_to_plot ) %>% 
      ggplot( aes(x = log10(weight_frac_observed) ,y = log10(weight_frac_predicted) ) ) + 
      geom_point(alpha=0.1) +
      geom_abline(slope = 1, intercept = 0 , col = "grey") +
      facet_grid( source~ mech_model  * growth_condition , scales="free" )+
      xlab("log10 observed metabolic proteome weight fraction")+
      ylab("log10 predicted metabolic proteome weight fraction") +
      geom_text(data = performance_df_s , 
                aes(x=pos_x, y=pos_y,
                    label = paste0( "n=",n,"\nR2=",round(R2,digits = 2), "\nRMSE=",round(RMSE,digits = 3)) ), 
                    size = 2.75 , 
                    alpha= 0.5 ) +
      theme(strip.text.y = element_text(size = 2))+
      theme_base() +
      theme(plot.background=element_rect( colour=NA))
    return(plot_out)
}