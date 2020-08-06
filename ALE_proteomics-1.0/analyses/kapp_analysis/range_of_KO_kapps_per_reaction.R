write_pdf <- FALSE

# investigate the spread of kapps across samples by comparing max to min



library(tidyverse)

load("output_data/kappmax/KO_MFA_kappmax/KO_MFA_kappmax.rdata")

t3_1tn_mfa <- res_li$t3_1tn_mfa

head( t3_1tn_mfa[ ! is.na(t3_1tn_mfa$kapp_per_s) ,  ] )

kapp_ranges <- t3_1tn_mfa %>% group_by( react_id ,sample_id) %>% 
  summarise( avg_kapp_overbrs = mean(kapp_per_s , na.rm = TRUE)) %>% # average kapps over bio reps
  group_by( react_id) %>% 
  mutate( n_kapp_per_react = sum(!is.na(avg_kapp_overbrs)) ) %>% 
  filter( n_kapp_per_react > 0 ) %>% 
  summarise( kapp_min = min(avg_kapp_overbrs , na.rm = TRUE) , 
             kapp_max = max(avg_kapp_overbrs , na.rm = TRUE) , 
             kapp_log10_sd = sd( log10(avg_kapp_overbrs) , na.rm = TRUE) , 
             n_kapp_per_react = first(n_kapp_per_react) ) %>% 
  select( react_id ,  kapp_min , kapp_max , n_kapp_per_react , kapp_log10_sd ) %>% 
  mutate( log10diff = log10(kapp_max/kapp_min) , 
          log2diff = log2(kapp_max/kapp_min) , 
          abs_diff = kapp_max - kapp_min ) %>% 
  arrange( desc(log10diff )) %>% 
  print(n=Inf)


kapp_var <- t3_1tn_mfa %>% group_by( react_id ,sample_id) %>% 
  summarise( avg_kapp_overbrs = mean(kapp_per_s , na.rm = TRUE)) %>% # average kapps over bio reps
  group_by( react_id) %>% 
  mutate( n_kapp_per_react = sum(!is.na(avg_kapp_overbrs)) ) %>% 
  filter( n_kapp_per_react > 0 ) %>% 
  summarise( kapp_min = min(avg_kapp_overbrs , na.rm = TRUE) , kapp_max = max(avg_kapp_overbrs , na.rm = TRUE) , n_kapp_per_react = first(n_kapp_per_react) )

# 
summary( kapp_ranges$log10diff[ kapp_ranges$n_kapp_per_react>10 ] )

summary( kapp_ranges$log2diff[ kapp_ranges$n_kapp_per_react>10 ] )

summary( kapp_ranges$kapp_log10_sd[ kapp_ranges$n_kapp_per_react>10 ] )

library(ggthemes)

kapp_hist <- kapp_ranges %>% filter(n_kapp_per_react>10) %>% 
  ggplot( aes( x = log2diff ) ) +
  geom_histogram()+
  #geom_density() +
  #stat_bin(geom="step") +
  xlim(0,NA) + 
  xlab( expression( "log"[2]*" "*frac("max("*italic(k)[app]*")","min("*italic(k)[app]*")")*" per reaction" ) ) +
  theme_base() +
  theme(plot.background=element_rect( colour=NA))

kapp_hist

if(write_pdf){
  ggsave(file = "Manuscript/Figures/kapp_structure/kapp_range_hist.pdf" , plot = kapp_hist , width = 7, height = 5.5 ,device = "pdf")
  kapp_hist
  dev.off()
}

###################################################################
# compare bilogical variation with variation due to KO+ALE

t3_1tn_mfa_s <- t3_1tn_mfa %>% select(react_id,sample_id,Master.Protein.Accessions,bio_rep , kapp_per_s)

t3_1tn_mfa_s_sj <- merge(t3_1tn_mfa_s , t3_1tn_mfa_s , by = c("react_id","sample_id") )

kapp_bio_rep_sd <- t3_1tn_mfa_s_sj %>% filter(bio_rep.x == "B1" & bio_rep.y == "B2") %>% 
  filter(! is.na(kapp_per_s.x), ! is.na(kapp_per_s.y) ) %>% 
  rowwise() %>% 
  mutate( #log2diff_reps = log2( max(kapp_per_s.x , kapp_per_s.y , na.rm = T) / min(kapp_per_s.x , kapp_per_s.y ,na.rm = T ) ) ,
    log10_sd = ( sd( c( log10(kapp_per_s.x) , log10(kapp_per_s.y) ) ) )
  ) %>% 
  #group_by(react_id) %>% 
  #summarise ( m.log2diff_reps = mean(log2diff_reps , na.rm = TRUE) ) %>% 
  ungroup() 

head(kapp_bio_rep_sd)
summary( kapp_bio_rep_sd$log10_sd )

summary( kapp_ranges$kapp_log10_sd )

wilc <- wilcox.test( kapp_bio_rep_sd$log10_sd , kapp_ranges$kapp_log10_sd )

length(kapp_bio_rep_sd$log10_sd[!is.na(kapp_bio_rep_sd$log10_sd)])
length(kapp_ranges$kapp_log10_sd[!is.na(kapp_ranges$kapp_log10_sd)])

wilc$statistic

