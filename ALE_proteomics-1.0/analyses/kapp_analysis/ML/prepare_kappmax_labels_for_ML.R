library(dplyr)

#############################################################################
#kappmax_KO_ALE
load("output_data/kappmax/KO_MFA_kappmax/KO_MFA_kappmax.rdata")
t3_1tn_mfa_kappmax <- as.data.frame(res_li$t3_1tn_mfa_kappmax)

# just keep reaction id and y
KO_mfa_kappmax <- subset(t3_1tn_mfa_kappmax , select = c(react_id , kappmax_per_s_br_avg) )

# give y name
names(KO_mfa_kappmax)[ names(KO_mfa_kappmax) == "kappmax_per_s_br_avg" ] <- "kappmax_KO_ALE_per_pp_per_s" 
# use standard name for react ids
names(KO_mfa_kappmax)[ names(KO_mfa_kappmax) == "react_id" ] <- "bigg_id" 
# remove NAs
KO_mfa_kappmax <- KO_mfa_kappmax[ ! is.na(KO_mfa_kappmax$kappmax_KO_ALE) , ]


saveRDS(KO_mfa_kappmax , file = "output_data/kappmax/KO_MFA_kappmax/KO_MFA_kappmax_Y.rds" ) 


#############################################################################
# kappmax Davidi
davidi <- read.csv("data/Davidi/Davidi_2016_S01_kmax.csv" , stringsAsFactors = FALSE, skip = 2)
head(davidi)

davidi$reaction..model.name.sybil <- sub( "(.+)_reverse$" ,"\\1_b", davidi$reaction..model.name. )

# check size of overlap:
length( setdiff( KO_mfa_kappmax$bigg_id , davidi$reaction..model.name.sybil))
length( intersect( KO_mfa_kappmax$bigg_id , davidi$reaction..model.name.sybil))
length( setdiff( davidi$reaction..model.name.sybil , KO_mfa_kappmax$bigg_id))

davidi_kappmax <- subset(davidi , select = c(reaction..model.name.sybil , kmax.per.polypeptide.chain..s.1.))

names(davidi_kappmax) [names(davidi_kappmax)  == "kmax.per.polypeptide.chain..s.1." ] <- "kappmax_davidi_per_pp_per_s" 
names(davidi_kappmax) [names(davidi_kappmax)  == "reaction..model.name.sybil" ] <- "bigg_id" 

head(davidi_kappmax)

saveRDS(davidi_kappmax , file = "output_data/kappmax/Davidi_kappmax/davidi_kappmax_Y.rds" ) 

#############################################################################
# maximum of KO ALE kappmax and Davidi et al kappmax

KO_mfa_kappmax_rbind <- KO_mfa_kappmax
davidi_kappmax_rbind <- davidi_kappmax 

names(KO_mfa_kappmax_rbind)[names(KO_mfa_kappmax_rbind) == "kappmax_KO_ALE_per_pp_per_s"] <- "kappmax_KO_ALE_davidi_per_pp_per_s"
names(davidi_kappmax_rbind)[names(davidi_kappmax_rbind) == "kappmax_davidi_per_pp_per_s"] <- "kappmax_KO_ALE_davidi_per_pp_per_s"

KO_mfa_davidi_kappmax <- rbind( KO_mfa_kappmax_rbind , 
                                davidi_kappmax_rbind ) 

# unique reactions with kappmax
length(unique(KO_mfa_davidi_kappmax$bigg_id))

# take the max per reaction
KO_mfa_davidi_kappmax <- KO_mfa_davidi_kappmax %>% 
  group_by(bigg_id) %>% 
  summarise( kappmax_KO_ALE_davidi_per_pp_per_s =  max(kappmax_KO_ALE_davidi_per_pp_per_s) )

stopifnot( ! any(is.na(KO_mfa_davidi_kappmax)) ) 

saveRDS(KO_mfa_davidi_kappmax , file = "output_data/kappmax/KO_MFA_Davidi_kappmax/KO_mfa_davidi_kappmax_Y.rds" ) 

#############################################################################
# in vitro kcat from kcat ML 

# features with kcat 
load_env <- new.env()

# features with kcat
load( "data/ML_inputs/ML_features_and_labels.rdata" ,
     envir = load_env
)
glob_kcat <- load_env$glob_feat_kcat_wo_exch

# bigg id s are duplicated because there are kcats for multiple assay conditions: summ,arise these 
glob_kcat_unique <- glob_kcat %>% group_by(bigg_id) %>% 
  summarise( log10_kcat_uni = ifelse( all(is.na(pH)) & all(is.na(one.Temp)) , 
                                      mean(log10_kcat) , 
                                      ifelse(all(is.na(one.Temp)),
                                             log10_kcat[which.min(pH_abs_dist_7.5)],
                                             log10_kcat[which.min(abs((1/one.Temp)-25))]
                                      )
  )
  )

kcat_in_vitro_ML <- glob_kcat_unique 

kcat_in_vitro_ML <- kcat_in_vitro_ML[ ! is.na(kcat_in_vitro_ML$log10_kcat_uni) , ]
kcat_in_vitro_ML$kcat_iv_ML_per_AS_per_s <- 10^kcat_in_vitro_ML$log10_kcat_uni
kcat_in_vitro_ML$bigg_id <- sub("_f$" , "", kcat_in_vitro_ML$bigg_id)

kcat_in_vitro_ML <-  subset ( kcat_in_vitro_ML , select = - c(log10_kcat_uni) )



head(kcat_in_vitro_ML)

saveRDS(kcat_in_vitro_ML , file = "output_data/kcat_in_vitro/kcat_in_vitro_ML_Y.rds" ) 

#############################################################################