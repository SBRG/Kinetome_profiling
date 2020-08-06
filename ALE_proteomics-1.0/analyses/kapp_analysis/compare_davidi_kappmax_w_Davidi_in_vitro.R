library(tidyverse)
library(ggrepel)
library(venn)
library(scales)

write_pdf <- FALSE

options(tibble.width = Inf,
        max.print = 1200
) # for printing tbls


davidi_in_vitro <- read.csv("data/Davidi/Davidi_2016_S01_kcat_vs_kmax.csv" , stringsAsFactors = FALSE, skip = 2)

davidi <- read.csv("data/Davidi/Davidi_2016_S01_kmax.csv" , stringsAsFactors = FALSE, skip = 2)

davidi$reaction..model.name.sybil <- sub( "(.+)_reverse$" ,"\\1_b", davidi$reaction..model.name. )
davidi_in_vitro$reaction..model.name.sybil <- sub( "(.+)_reverse$" ,"\\1_b", davidi_in_vitro$reaction..model.name. )

#####################################################


davidi_viv_vit <- merge( davidi , davidi_in_vitro , by = "reaction..model.name.sybil" )


##############################################################

cortest <- with( davidi_viv_vit , cor.test( log10(kmax.per.polypeptide.chain..s.1.) , log10(kcat..s.1. )  )  )
cortest$estimate^2

MAE <- with( davidi_viv_vit , mean( abs( log10(kmax.per.polypeptide.chain..s.1. ) - log10(kcat..s.1.) ) )  )

# count complete obs 
n_compl <- davidi_viv_vit %>% 
  summarise( n_compl =  sum(! is.na(log10(kmax.per.polypeptide.chain..s.1. ) & ! is.na( log10(kcat..s.1.) ) )) ) %>% unlist()

# scatterplot
ggp <- davidi_viv_vit %>% ggplot( aes(x = log10(kmax.per.polypeptide.chain..s.1. ) , log10(kcat..s.1.) )) + 
  geom_abline(slope = 1 , intercept= c(-1,0,1) , col = c("grey","black","grey") ) +
  #geom_text_repel( aes(label = react_id) , size = 2.6, alpha = 0.3, segment.alpha = 0.1) +
  geom_text_repel( aes(label = reaction..model.name.sybil ) , 
                   size = 2.6, alpha = 0.2, 
                   segment.alpha = 0.1 ,  
                   data = davidi_viv_vit[ abs( log10(davidi_viv_vit$kmax.per.polypeptide.chain..s.1.) - log10(davidi_viv_vit$kcat..s.1.) ) > 1 ,  ] ,   
                   inherit.aes = T) +
  geom_point( alpha =0.3) + ggthemes::theme_base() + 
  theme(panel.grid = element_line(color = "lightgrey",linetype =3) , 
        plot.background=element_rect(fill="white", colour=NA) , 
        #legend.background = element_blank(),
        legend.background = element_rect(fill=alpha('white', 0.1),colour = "black"),
        #legend.box.background = element_rect(colour = "black"),
        legend.position = c(0.8, 0.15),
        legend.title = element_blank()
  )  +
  #xlab("log10( kappmax in Davidi et al. )") + 
  #ylab("log10( kcat in vitro )") +
  xlab(expression("log"[10]*"( "*italic("k")["app,max"]*" from growth conditions )") ) +
  ylab(expression("log"[10]*"( "*italic("k")["cat"]*" in vitro )") ) +
  annotate("text" ,x = -2 , y =2.2 , parse = TRUE , size = 5 , 
           label = paste0( "italic(R)^2 ==" ,round( cortest$estimate^2 , 2) )
  ) +
  annotate("text" ,x = -2 , y =1.8 , parse = TRUE , size = 5 , 
           label = paste0( "MAE ==" ,round( MAE , 2) )
  ) +
  annotate("text" ,x = -2 , y =1.45 , parse = TRUE , size = 5 , 
           label = paste0("n == ", n_compl ) )


ggp



if(write_pdf){
  ggsave(file = "Manuscript/Figures/kappmax_comparison/kappmax_davidi_vs_kcat_davidi.pdf" , plot = ggp , width = 5, height = 5 ,device = "pdf")
}
