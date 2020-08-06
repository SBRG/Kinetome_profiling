

barplot_mechmod_perf <- function(performance_df , # source , mech_model, condition, R2, RMSE
                                 mods_to_to_plot, # which models with keff sources listed in "source"
                                 rename_map, # 
                                 fill_cols, legend_labels , legend_title
                                 ) {
  
    library(ggthemes)
    library(ggplot2)
    library(gridExtra)
    library(scales)
    
    library(dplyr)
    
    table(performance_df$source)
    
    performance_df %>% group_by(source , mech_model) %>% summarise( m.R2 = mean(R2) , m.RMSE = mean(RMSE))
    

    performance_df_filt <- performance_df %>% 
      filter(source %in% mods_to_to_plot )
    

    performance_df_filt$source <- factor(performance_df_filt$source , levels = mods_to_to_plot )
    
    performance_df_filt$condition <- paste0(performance_df_filt$condition, "\n(n=" , performance_df_filt$n , ")" )
    

    ##############
    
    barp_lwd <- 0.5
    
    R2_MOM <- performance_df_filt %>% 
      filter(mech_model=="MOMENT") %>%
      ggplot( aes( x = condition , y = R2 , fill = source) ) + 
      geom_bar(stat="identity", position=position_dodge() , lwd = barp_lwd ) +
      ggthemes::theme_few() + 
      ggthemes::scale_colour_few()+
      #facet_grid(.~mech_model) +
      ylab (expression("R "^2) )
    
    RMSE_MOM <- performance_df_filt %>% 
      filter(mech_model=="MOMENT") %>%
      ggplot( aes( x = condition , y = RMSE , fill = source) ) + 
      geom_bar(stat="identity", position=position_dodge(), lwd = barp_lwd)+
      ggthemes::theme_few() + 
      #facet_grid(.~mech_model) +
      ggthemes::scale_colour_few() +
      theme( panel.border = element_rect( size=0.5) )
    
    #############
    
    
    R2_ME <- performance_df_filt %>% 
      filter(mech_model=="ME") %>%
      #select( - mech_model) %>%
      ggplot( aes( x = condition , y = R2 , fill = source) ) + 
      geom_bar(stat="identity", position=position_dodge(), lwd = barp_lwd) +
      ggthemes::theme_few() + 
      ggthemes::scale_colour_few()+
      #facet_grid(.~mech_model) +
      ylab (expression("R "^2) )
    
    RMSE_ME <- performance_df_filt %>% 
      filter(mech_model=="ME") %>%
      #select( - mech_model) %>%
      ggplot( aes( x = condition , y = RMSE , fill = source) ) + 
      geom_bar(stat="identity", position=position_dodge(), lwd = barp_lwd)+
      ggthemes::theme_few() + 
      #facet_grid(.~mech_model) +
      ggthemes::scale_colour_few()  +
      theme( panel.border = element_rect( size=0.5) )

    
    ##################################################################################################################
    # merge plots 
    
    # extract single legend 
    g <- ggplotGrob( R2_MOM + 
                       theme(legend.position="bottom", legend.title.align=0.5 , text = element_text(size =9) , legend.text.align = 0 , 
                             legend.key.size = unit(3,"mm") ) + 
                       scale_fill_manual(values=fill_cols , 
                                          name = ,legend_title,
                                          labels = legend_labels
                                         )
                       #) #+ guides(fill=guide_legend(nrow=2,byrow=TRUE))
                     #scale_fill_discrete()
    )$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]] 

    arranged_grob <- grid.arrange( arrangeGrob( RMSE_MOM + 
                                 theme(legend.position="none" ,  axis.title.x=element_blank() , plot.title = element_text(hjust = 0.5),text = element_text(size =9) ) + 
                                 ggtitle("MOMENT model predictions") + 
                                 scale_fill_manual(values=fill_cols) +
                                 ylab("RMSE (on log10 weight fractions)")
                               ,
                               RMSE_ME + 
                                 theme(legend.position="none", #axis.title.y=element_blank(), 
                                       plot.title = element_text(hjust = 0.5) , text = element_text(size =9) #, axis.title.x=element_blank() 
                                 ) +
                                 ggtitle("ME model predictions") + 
                                 xlab("growth substrate\n(number of comparisons)") + 
                                 ylab("RMSE (on log10 weight fractions)") +
                                 scale_fill_manual(values=fill_cols)
    ),
    legend ,
    ncol = 1,
    heights = c(8,1)
    )
    
    arranged_grob
    
    
    return( list ( arranged_grob = arranged_grob) ) 
}