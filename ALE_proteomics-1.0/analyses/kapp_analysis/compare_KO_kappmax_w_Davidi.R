library(tidyverse)
library(ggrepel)
library(venn)
library(scales)

write_pdf <- FALSE

options(tibble.width = Inf,
        max.print = 1200
) # for printing tbls

# created in "calc_kapp_MFA.R"
t3_1tn_mfa <- read.csv("output_data/kappmax/KO_MFA_kappmax/KO_MFA_kappmax_full.csv" , stringsAsFactors = FALSE )
t3_1tn_mfa_kappmax <- read.csv("output_data/kappmax/KO_MFA_kappmax/KO_MFA_kappmax_summarised.csv" , stringsAsFactors = FALSE )

t3_1tn_mfa_kappmax <- t3_1tn_mfa_kappmax[! is.na(t3_1tn_mfa_kappmax$kappmax_per_s_br_avg) , ]

davidi_in_vitro <- read.csv("data/Davidi//Davidi_2016_S01_kcat_vs_kmax.csv" , stringsAsFactors = FALSE, skip = 2)

davidi <- read.csv("data/Davidi/Davidi_2016_S01_kmax.csv" , stringsAsFactors = FALSE, skip = 2)
head(davidi)

davidi$reaction..model.name.sybil <- sub( "(.+)_reverse$" ,"\\1_b", davidi$reaction..model.name. )
davidi_in_vitro$reaction..model.name.sybil <- sub( "(.+)_reverse$" ,"\\1_b", davidi_in_vitro$reaction..model.name. )

#####################################################

# check size of overlap:
length( setdiff( t3_1tn_mfa_kappmax$react_id , davidi$reaction..model.name.sybil))
length( intersect( t3_1tn_mfa_kappmax$react_id , davidi$reaction..model.name.sybil))
length( setdiff( davidi$reaction..model.name.sybil , t3_1tn_mfa_kappmax$react_id))
# plot venn to illustrate overlap
venn( list(t3_1tn_mfa_kappmax$react_id , davidi$reaction..model.name.sybil) , zcolor = "style" , snames = c( NA  , NA), cexil = 2.5   )

if(write_pdf){
  pdf(file = "Manuscript/Figures/kappmax_comparison/kappmax_KO_vs_kappmax_davidi_venn.pdf" , width = 5, height = 5 )
  venn( list(t3_1tn_mfa_kappmax$react_id , davidi$reaction..model.name.sybil) , zcolor = "style" , snames = c( NA  , NA)  , cexil = 2.5 )
  dev.off()
}


t3_1tn_mfa_kappmax_dav <- merge( t3_1tn_mfa_kappmax , davidi , by.x = "react_id" ,by.y =  "reaction..model.name.sybil" )

t3_1tn_mfa_kappmax_dav <- t3_1tn_mfa_kappmax_dav %>% 
  mutate( sig_diff = kmax.per.polypeptide.chain..s.1. < kappmax_samp_95CI_lb | kmax.per.polypeptide.chain..s.1. > kappmax_samp_95CI_ub )

t3_1tn_mfa_kappmax_dav %>% filter(sig_diff) %>% mutate( sign = ifelse(kappmax_per_s_br_avg <kmax.per.polypeptide.chain..s.1. , "smaller","greater")) %>% 
  group_by(sign) %>%  summarise( n = n())

# count non-significant cases for main text 
table(t3_1tn_mfa_kappmax_dav$sig_diff)

##############################################################
# correlation of kappmax in KO ALE vs Davidi

cortest <- with( t3_1tn_mfa_kappmax_dav , cor.test(log10(kappmax_per_s_br_avg ) , log10(kmax.per.polypeptide.chain..s.1.) )  )
cortest$estimate^2

MAE <- with( t3_1tn_mfa_kappmax_dav , mean( abs( log10(kappmax_per_s_br_avg ) - log10(kmax.per.polypeptide.chain..s.1.) ) )  )


# corr test for non-significant differences
sum(!t3_1tn_mfa_kappmax_dav$sig_diff)
# 66 non-sig
cortest_ndiff <- with( t3_1tn_mfa_kappmax_dav[ !t3_1tn_mfa_kappmax_dav$sig_diff , ] , cor.test(log10(kappmax_per_s_br_avg ) , log10(kmax.per.polypeptide.chain..s.1.) )  )
cortest_ndiff$estimate^2

# count complete obs 
n_compl <- t3_1tn_mfa_kappmax_dav %>% 
  summarise( n_compl =  sum(! is.na(log10(kappmax_per_s_br_avg ) & ! is.na( log10(kmax.per.polypeptide.chain..s.1.) ) )) ) %>% unlist()

# scatterplot
ggp <- t3_1tn_mfa_kappmax_dav %>% ggplot( aes(x = log10(kappmax_per_s_br_avg ) , log10(kmax.per.polypeptide.chain..s.1.) )) + 
  geom_abline(slope = 1 , intercept= c(-1,0,1) , col = c("grey","black","grey") ) +
  geom_errorbarh(aes(xmin = log10( kappmax_samp_95CI_lb) , 
                     xmax = log10( kappmax_samp_95CI_ub) ),
                 alpha =0.15 , size = 0.4) +
  #geom_text_repel( aes(label = react_id) , size = 2.6, alpha = 0.3, segment.alpha = 0.1) +
  geom_text_repel( aes(label = react_id) , 
                   size = 2.6, alpha = 0.2, 
                   segment.alpha = 0.1 ,  
                   data = t3_1tn_mfa_kappmax_dav[ abs( log10(t3_1tn_mfa_kappmax_dav$kappmax_per_s_br_avg) - log10(t3_1tn_mfa_kappmax_dav$kmax.per.polypeptide.chain..s.1.) ) > 1 ,  ] ,   
                   inherit.aes = T) +
  geom_point(aes( col = sig_diff ), alpha =0.6) + ggthemes::theme_base() + 
  theme(panel.grid = element_line(color = "lightgrey",linetype =3) , 
        plot.background=element_rect(fill="white", colour=NA) , 
        #legend.background = element_blank(),
        legend.background = element_rect(fill=alpha('white', 0.1),colour = "black"),
        #legend.box.background = element_rect(colour = "black"),
        legend.position = c(0.8, 0.15),
        legend.title = element_blank()
        )  +
  xlab(expression("log"[10]*"( "*italic("k")["app,max"]*" from KO ALEs )") ) +
  ylab(expression("log"[10]*"( "*italic("k")["app,max"]*" from growth conditions )") ) +
  scale_color_discrete(name = "", labels = c("non-significant", "significant")) +
  scale_x_continuous(breaks = pretty(log10(t3_1tn_mfa_kappmax_dav$kappmax_per_s_br_avg ), n = 5)) +
  scale_y_continuous(breaks = pretty(log10(t3_1tn_mfa_kappmax_dav$kmax.per.polypeptide.chain..s.1.), n = 5)) +
  annotate("text" ,x = -3 , y =2.2 , parse = TRUE , size = 5 , 
           label = paste0( "italic(R)^2 ==" ,round( cortest$estimate^2 , 2) )
  ) +
  annotate("text" ,x = -3 , y =1.7 , parse = TRUE , size = 5 , 
           label = paste0( "MAE ==" ,round( MAE , 2) )
  ) +
  annotate("text" ,x = -3 , y =1.2 , parse = TRUE , size = 5 , 
           label = paste0("n == ", n_compl ) )




ggp


if(write_pdf){
  ggsave(file = "Manuscript/Figures/kappmax_comparison/kappmax_KO_vs_kappmax_davidi.pdf" , plot = ggp , width = 5, height = 5 ,device = "pdf")
}
