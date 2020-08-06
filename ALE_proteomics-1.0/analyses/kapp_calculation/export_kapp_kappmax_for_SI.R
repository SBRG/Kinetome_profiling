library(dplyr)

######################################################

t3_1tn_mfa <- read.csv("output_data/kappmax/KO_MFA_kappmax/KO_MFA_kappmax_full.csv" , stringsAsFactors = FALSE)

head(t3_1tn_mfa[ ! is.na(t3_1tn_mfa$mfa_sampling_ave ),])

t3_1tn_mfa_SI <- t3_1tn_mfa %>% select(react_id , direction, sample_id,Master.Protein.Accessions,bnum, bio_rep , fmol_per_gDW_calib_avg,
                      flux_units, mfa_sampling_ave, mfa_sampling_var, mfa_sampling_lb, mfa_sampling_ub, # MFA data
                      kapp_per_s #, kappmax_per_s_br, kappmax_per_s_br_n , kappmax_per_s_br_sample_id_of_max ,kappmax_per_s_br_avg , kappmax_per_s_br_sd
                      )

t3_1tn_mfa_SI$sample_id <- sub("^EVO(1|2)$", "WT\\1" , t3_1tn_mfa_SI$sample_id)
t3_1tn_mfa_SI$sample_id <- sub("EVO|Evo", "" , t3_1tn_mfa_SI$sample_id)

head(t3_1tn_mfa_SI)
tail(t3_1tn_mfa_SI)

# react_id : BIGG reaction id 
# direction : Reaction direction relative to BIGG definition
# sample_id : ID of ALE endpoint strain
# Master.Protein.Accessions : gene name
# bnum : gene "b-number"
# bio_rep :  Biological replicate number
# fmol_per_gDW_calib_avg : Protein abundance in fmol gDW-1
# flux_units : Unit of flux estimate
# mfa_sampling_ave : MFA flux sampling average
# mfa_sampling_var : MFA flux sampling variance
# mfa_sampling_lb : MFA flux lower bound
# mfa_sampling_ub : MFA flux upper bound
# kapp_per_s : kapp in s-1

write.csv( t3_1tn_mfa_SI , file = "output_data/SM/protein_ab_flux_kapp_table/protein_abundance_flux_kapp_table.csv" , row.names = FALSE)

######################################################

t3_1tn_mfa_kappmax <-  read.csv( "output_data/kappmax/KO_MFA_kappmax/KO_MFA_kappmax_summarised.csv" , stringsAsFactors = FALSE )

head(t3_1tn_mfa_kappmax)
t3_1tn_mfa_kappmax_SI <-  t3_1tn_mfa_kappmax %>% select( - X , -all_r_per_g_homo , -one_to_one) %>%  filter(! is.na(kappmax_per_s_br_avg) )

t3_1tn_mfa_kappmax_SI$kappmax_per_s_br_avg_sample_id_of_max <- sub("^EVO(1|2)$", "WT\\1" , t3_1tn_mfa_kappmax_SI$kappmax_per_s_br_avg_sample_id_of_max)
t3_1tn_mfa_kappmax_SI$kappmax_per_s_br_avg_sample_id_of_max <- sub("EVO|Evo", "" , t3_1tn_mfa_kappmax_SI$kappmax_per_s_br_avg_sample_id_of_max)

t3_1tn_mfa_kappmax_SI$kappmax_samp_maj_sample_id <- sub("^EVO(1|2)$", "WT\\1" , t3_1tn_mfa_kappmax_SI$kappmax_samp_maj_sample_id)
t3_1tn_mfa_kappmax_SI$kappmax_samp_maj_sample_id <- sub("EVO|Evo", "" , t3_1tn_mfa_kappmax_SI$kappmax_samp_maj_sample_id)

head(t3_1tn_mfa_kappmax_SI)

# react_id : BIGG reaction ID, where _b indicates the opposite reaction direction
# Master.Protein.Accessions : Gene name
# kappmax_per_s_br_avg : kapp,max, averaged over biological replicates
# kappmax_per_s_br_sd : Standard deviation of kappmax over biological replicates
# kappmax_per_s_br_avg_n : Number of biological replicates available
# kappmax_per_s_n_avg_brs_used : Number of samples used to calculate kapp,max
# kappmax_per_s_br_avg_sample_id_of_max : Sample id that contributed the highest kapp
# kappmax_samp_avg : Average kapp,max in parametric bootstrap
# kappmax_samp_med : Median kapp,max in parametric bootstrap
# kappmax_samp_sd : Standard deviation of kapp,max in parametric bootstrap
# kappmax_samp_95CI_lb : Lower 95% confidence interval of kapp,max in parametric bootstrap
# kappmax_samp_95CI_ub : Upper 95% confidence interval of kapp,max in parametric bootstrap
# kappmax_samp_maj_sample_id : Sample id that contributed the highest kapp in parametric bootstrap

write.csv( t3_1tn_mfa_kappmax_SI , file = "output_data/SM/kappmax_table/kappmax_table.csv" , row.names = FALSE)

