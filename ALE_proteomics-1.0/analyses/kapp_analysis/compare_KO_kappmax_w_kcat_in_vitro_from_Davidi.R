library(tidyverse)
library("ggrepel")

write_pdf <- FALSE

options(tibble.width = Inf,
        max.print = 1200
) # for printing tbls

# created in "calc_kapp_MFA.R"
t3_1tn_mfa <- read.csv("output_data/kappmax/KO_MFA_kappmax/KO_MFA_kappmax_full.csv" , stringsAsFactors = FALSE )
t3_1tn_mfa_kappmax <- read.csv("output_data/kappmax/KO_MFA_kappmax/KO_MFA_kappmax_summarised.csv" , stringsAsFactors = FALSE )
t3_1tn_mfa_kappmax <- t3_1tn_mfa_kappmax[! is.na(t3_1tn_mfa_kappmax$kappmax_per_s_br_avg) , ]
# use only 1t1 cases:
#t3_1tn_mfa_kappmax <- t3_1tn_mfa_kappmax[ t3_1tn_mfa_kappmax$one_to_one , ]

davidi_in_vitro <- read.csv("data/Davidi/Davidi_2016_S01_kcat_vs_kmax.csv" , stringsAsFactors = FALSE, skip = 2)
davidi_in_vitro$reaction..model.name.sybil <- sub( "(.+)_reverse$" ,"\\1_b", davidi_in_vitro$reaction..model.name. )

#############################

length(intersect(davidi_in_vitro$reaction..model.name.sybil , t3_1tn_mfa_kappmax$react_id ))

# inner join
t3_1tn_mfa_kappmax_davidi_vitro <- merge( t3_1tn_mfa_kappmax , davidi_in_vitro , by.x = "react_id" , by.y = "reaction..model.name.sybil")

###############################

# add significant differences based on confidence intervals
t3_1tn_mfa_kappmax_davidi_vitro <- t3_1tn_mfa_kappmax_davidi_vitro %>% mutate( sig_diff = kcat..s.1. < kappmax_samp_95CI_lb | kcat..s.1. > kappmax_samp_95CI_ub )

cortest <- with( t3_1tn_mfa_kappmax_davidi_vitro , cor.test(log10(kappmax_per_s_br_avg ) , log10(kcat..s.1.) )  )
cortest$estimate^2

MAE <- with( t3_1tn_mfa_kappmax_davidi_vitro , mean( abs( log10(kappmax_per_s_br_avg ) - log10(kcat..s.1.) ) )  )

# count greater and smaller cases
t3_1tn_mfa_kappmax_davidi_vitro %>% filter(sig_diff) %>% mutate( sign = ifelse(kappmax_per_s_br_avg < kcat..s.1. , "smaller","greater")) %>% 
  group_by(sign) %>%  summarise( n = n())

# fraction of in vitro kcats that fall into confidence intervals
t3_1tn_mfa_kappmax_davidi_vitro %>% summarize(sum(!sig_diff) / n() )
#  25% of in vitro kcats fall into confidence intervals. Goelzer 2005 found 65%, but their CIs are much wider as well, so these results are not comparable

# count complete obs 
n_compl <- t3_1tn_mfa_kappmax_davidi_vitro %>% summarise( n_compl =  sum(! is.na(log10(kcat..s.1.) & ! is.na( log10(kappmax_per_s_br_avg ) ) )) )

# count non-significant cases for main text 
table(t3_1tn_mfa_kappmax_davidi_vitro$sig_diff)

##########################

# scatterplot
ggp <- t3_1tn_mfa_kappmax_davidi_vitro %>% ggplot( aes( x =log10(kappmax_per_s_br_avg )  , y = log10(kcat..s.1.)  )) + 
  geom_abline(slope = 1 , intercept= c(-1,0,1) , col = c("grey","black","grey") ) +
  geom_errorbarh(aes(xmin = log10( kappmax_samp_95CI_lb) , 
                    xmax = log10( kappmax_samp_95CI_ub) ),
                alpha =0.15)+
  #geom_text_repel( aes(label = gene.name) , size = 3, alpha = 0.3, segment.alpha = 0.1) +
  geom_point( aes(col = sig_diff ) , alpha =0.6 )+ 
  geom_text_repel( aes(label = react_id) , 
                   size = 2.6, alpha = 0.2, 
                   segment.alpha = 0.1 ,  
                   data = t3_1tn_mfa_kappmax_davidi_vitro[ abs( log10(t3_1tn_mfa_kappmax_davidi_vitro$kappmax_per_s_br_avg) - log10(t3_1tn_mfa_kappmax_davidi_vitro$kcat..s.1.) ) > 1 ,  ] ,   
                   inherit.aes = T)+ ggthemes::theme_base() + 
  theme(panel.grid = element_line(color = "lightgrey",linetype =3) , 
        plot.background=element_rect(fill="white", colour=NA) , 
        #legend.background = element_blank(),
        legend.background = element_rect(fill=alpha('white', 0.1),colour = "black"),
        #legend.box.background = element_rect(colour = "black"),
        legend.position = c(0.8, 0.15),
        legend.title = element_blank()
  )  +
  scale_color_discrete(name = "", labels = c("non-significant", "significant")) +
  xlab(expression("log"[10]*"( "*italic("k")["app,max"]*" from KO ALEs )") ) +
  ylab(expression("log"[10]*"( "*italic("k")["cat"]*" in vitro )") ) +
  annotate("text" ,x = -2 , y =2.2 , parse = TRUE , size = 5 , 
           label = paste0( "italic(R)^2 ==" ,round( cortest$estimate^2 , 2) )
  ) +
  annotate("text" ,x = -2 , y =1.8 , parse = TRUE , size = 5 , 
           label = paste0( "MAE ==" ,round( MAE , 2) )
  ) +
  annotate("text" ,x = -2 , y =1.45 , parse = TRUE , size = 5 , 
           label = paste0("n == ", n_compl ) )
#ggthemes::theme_base() + theme(panel.grid = element_line(color = "lightgrey",linetype =3))


ggp

if(write_pdf){
  ggsave(file = "Manuscript/Figures/kappmax_comparison/kappmax_KO_vs_kcat_davidi.pdf" , plot = ggp , width = 5, height = 5 ,device = "pdf")
}

