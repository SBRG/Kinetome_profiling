# Parameters:
# fixed gP per gDW parameter
const_gP_per_gDW <- 0.55
# from Neidhardt:
# (2.35e6 molecs cell-1 / 6.0221e23 molecs per mol )* 1e15 fmol / mol  * 1/(2.8e-13 gDW/cell)
const_fmol_per_gDW <- (2.35e6/6.0221e23)*1e15 * 1/(2.8e-13)

library(MASS)
library(dplyr)
library(ggplot2)
library(broom)

# created by "parse_PSMs.R"
load("data/proteomics/ALE_Proteomics_1_26_2019.rdata")
ups2 <- read.csv("data/proteomics/UPS2.csv",stringsAsFactors = FALSE)
# the Amount.fmoles are for  1 vial of UPS2, which contains 10.6ug protein, of which 0.86ug were loaded per sample
ups2$UniProt_Accession_ext <- paste0( ups2$UniProt_Accession , "ups" )
 # mols * 10e-15 * g * 6.0221e23 * 
ups2$ups.weight_ug <- ups2$Amount_fmoles * 1e-15 *  ups2$MW_Da * 1e6 
sum(ups2$ups.weight_ug)
# 10.69ug, matches UPS manual (10.6)
# 0.88ug were added to 6ug sample, of which 1ug was finally loaded:
# amount in one batch of standard * fraction of batch spiked * fraction of spiked sample loaded 
ups2$ups.weight_loaded_ug <- ups2$ups.weight_ug * (0.88/sum(ups2$ups.weight_ug)) * 1/6.88
ups2$ups_Amount_fmoles_loaded <- ups2$Amount_fmoles * (0.88/sum(ups2$ups.weight_ug) ) * 1/6.88

# number of theoretically observble peptides created by "analyses/auxilliary_scripts/calc_n_observable_peptides.R"
n_obs_peps <- read.csv( file = "data/ALEdb_genomes/variant_files/IND_MG1655_F001/n_observable_peptides.csv" , stringsAsFactors = FALSE ) 

########################################################
# top3

# Silva et al about top3:
# "The average MS signal response for the three most intense tryptic peptides is calculated for each well characterized  protein in the mixture, 
# including those to the internal standard protein(s)."
# Fabre et al "The TOP3 is calculated as the mean of the three highest peptides areas measured for each protein"

top3_by_tech_r <- full_df_s %>% group_by( sample_id , bio_rep , Annotated.Sequence , Master.Protein.Accessions , Spectrum.File ) %>% # master or prot ?. Master makes more sense to me 
  summarise( sum_Area = sum ( Area , na.rm = TRUE ) ) %>% # sum up areas per peptide per protein per experiment
  group_by( sample_id  , bio_rep ,  Master.Protein.Accessions , Spectrum.File ) %>% # master or prot ?
  top_n(n = 3 , wt = sum_Area ) %>% 
  summarise( top3_tr =mean(sum_Area , na.rm = TRUE) ) %>% ungroup() %>%
  as.data.frame()

# calc average of top3 across techical replicates, and sd on log10 scale 
top3 <- top3_by_tech_r %>% group_by( sample_id , bio_rep , Master.Protein.Accessions ) %>%
  summarise( top3  = mean(top3_tr , na.rm = TRUE) , 
             top3_tech_sd_log10  = sd( log10(top3_tr) , na.rm = TRUE) )

ab_proxies <- top3 

############################################################
# calibration regression
# regress upsumol loaded on the top3
# fit individual models

top3_by_tech_r_ups <- merge( top3_by_tech_r , ups2, by.x = "Master.Protein.Accessions" , by.y = "UniProt_Accession_ext" ) 

# fit a OLS and a robust lm for log10(fmol abundance)
ups_lm_models <- top3_by_tech_r_ups %>% group_by(sample_id,bio_rep) %>% filter(top3_tr>0) %>%
     do( lm_model = lm( log10(ups_Amount_fmoles_loaded) ~ log10(top3_tr) , data = . ) , 
         rlm_model = rlm( log10(ups_Amount_fmoles_loaded) ~ log10(top3_tr) , data = . ) )

# inspect models
tidy( ups_lm_models , lm_model )
tidy( ups_lm_models , rlm_model )

# predict the ups 
ups_lm_preds <- ups_lm_models %>% augment(lm_model)
ups_rlm_preds <- ups_lm_models %>% augment(rlm_model)

# plot predictions for each model
# first, lm
ups_lm_preds %>% ggplot(aes( x = log10.top3_tr. , y=log10.ups_Amount_fmoles_loaded. )) + 
  geom_point(alpha = 0.2) + geom_line(aes(y = .fitted)) +
  facet_wrap( . ~ sample_id+bio_rep , ncol = 13)
# then, rlm
ups_rlm_preds %>% ggplot(aes( x = log10.top3_tr. , y=log10.ups_Amount_fmoles_loaded. )) + 
  geom_point(alpha = 0.2) + geom_line(aes(y = .fitted)) +
  facet_wrap( . ~ sample_id+bio_rep , ncol = 13)

#################################
# add ups models to abundance proxies and predict

top3_by_tech_r_w_mods <- as_tibble(top3_by_tech_r) %>% left_join( ups_lm_models , by = c("sample_id", "bio_rep" ) )
# add calbrated abundacne predictions to top3 data
top3_by_tech_r_w_pred <- top3_by_tech_r_w_mods %>% group_by(sample_id,bio_rep) %>% 
  do(modelr::add_predictions(., first(.$rlm_model) , var = "log10_fmol_pred")) %>% ungroup() %>%
  mutate( fmol_pred = 10^log10_fmol_pred ) 

####
# sanity check (calculate predictions by hand):
# first row
top3_by_tech_r_w_pred$rlm_model[[1]]$coefficients[1] + top3_by_tech_r_w_pred$rlm_model[[1]]$coefficients[2] * log10(top3_by_tech_r_w_pred$top3_tr[1])
# last row
top3_by_tech_r_w_pred$rlm_model[[nrow(top3_by_tech_r_w_pred)]]$coefficients[1] + top3_by_tech_r_w_pred$rlm_model[[nrow(top3_by_tech_r_w_pred)]]$coefficients[2] * log10(top3_by_tech_r_w_pred$top3_tr[nrow(top3_by_tech_r_w_pred)]) 
tail(top3_by_tech_r_w_pred)
# ok

# remove models form prediction dataframe foe easier handling:
top3_by_tech_r_w_pred <- top3_by_tech_r_w_pred %>% select(- lm_model, - rlm_model )

# plot abundacne dists
top3_by_tech_r_w_pred %>% ggplot(aes(x=log10_fmol_pred)) + geom_histogram() + facet_wrap(. ~sample_id + bio_rep)

#########################################
# convert units of predictions

# convert fmol abundaces to masses:
Glu4_masses <- read.csv("data/ALEdb_genomes/variant_files/GLU_A4_-_ref_or_evo04/GLU_A4_protein_masses.csv" , stringsAsFactors = FALSE)
setdiff( top3_by_tech_r_w_pred$Master.Protein.Accessions , Glu4_masses$gene_name )
# only UPS masses missing
# add masses to predictions
top3_by_tech_r_w_pred <- merge(top3_by_tech_r_w_pred  , Glu4_masses , by.x = "Master.Protein.Accessions" , by.y = "gene_name" , all.x = TRUE )
# convert fol to masses
top3_by_tech_r_w_pred$ug_pred <- top3_by_tech_r_w_pred$fmol_pred * 1e-15 * top3_by_tech_r_w_pred$mass_Da * 1e6 

# remove ups 
top3_by_tech_r_w_pred <- top3_by_tech_r_w_pred[ ! grepl("ups" , top3_by_tech_r_w_pred$Master.Protein.Accessions) , ]


# compute mass and amount sums per tech rep  
top3_by_tech_r_w_pred <- top3_by_tech_r_w_pred %>% group_by(Spectrum.File , sample_id , bio_rep) %>% 
  mutate(tech_rep_sum_fmol_pred = sum( fmol_pred ) ,
         tech_rep_frac_fmol_pred = fmol_pred / sum( fmol_pred ),
         tech_rep_sum_ug_pred = sum(ug_pred),
         tech_rep_frac_ug_pred = ug_pred/ sum(ug_pred),
         tech_rep_gP_per_gDW_rel = (ug_pred/ sum(ug_pred)) * const_gP_per_gDW ,
         tech_rep_fmol_per_gDW_calib = (fmol_pred / sum( fmol_pred )) * const_fmol_per_gDW,  # total protein from calibrated abundnces
         tech_rep_ug_per_gDW_calib = (ug_pred / sum( ug_pred )) * const_gP_per_gDW # total protein from calibrated abundnces
         ) 
###########################
# average predicted amounts across tech reps (and sd)
top3_preds <- top3_by_tech_r_w_pred %>% filter(top3_tr>0) %>% group_by( sample_id , bio_rep , Master.Protein.Accessions ) %>%
  summarise( fmol_pred_avg  = mean(fmol_pred , na.rm = TRUE) , 
             fmol_pred_avg_log10  = mean(log10_fmol_pred, na.rm = TRUE) , 
             fmol_pred_techsd_log10  = sd( log10(fmol_pred) , na.rm = TRUE),
             fmol_pred_techsd  = sd( fmol_pred , na.rm = TRUE),
             ug_pred_avg  = mean(ug_pred , na.rm = TRUE) , 
             ug_pred_med  = median(ug_pred , na.rm = TRUE) , 
             ug_pred_techsd_log10  = sd( log10(ug_pred) , na.rm = TRUE),
             ug_pred_techsd  = sd( ug_pred , na.rm = TRUE),
             fmol_per_gDW_calib_avg = mean( tech_rep_fmol_per_gDW_calib , na.rm=TRUE ),
             fmol_per_gDW_calib_techsd_log10 = sd( log10(tech_rep_fmol_per_gDW_calib) , na.rm=TRUE ),
             fmol_per_gDW_calib_techsd = sd( tech_rep_fmol_per_gDW_calib , na.rm=TRUE ) , 
             ug_per_gDW_calib_avg = mean(tech_rep_ug_per_gDW_calib, na.rm = TRUE),
             ug_per_gDW_calib_techsd_log10 = sd(log10(tech_rep_ug_per_gDW_calib), na.rm = TRUE),
             ug_per_gDW_calib_techsd = mean(tech_rep_ug_per_gDW_calib, na.rm = TRUE)
             )
# NAs in sd possible when  protein was detected in <2 reps

###############
# check sums of abundances per sample
top3_preds_sums <- top3_preds %>% group_by( sample_id , bio_rep) %>% 
  summarise(sum_fmol_pred_avg = sum(10^fmol_pred_avg_log10) ,
            sum_ug_pred_avg = sum(ug_pred_avg,na.rm = TRUE),
            sum_ug_pred_med = sum(ug_pred_med,na.rm = TRUE)
            )

top3_preds <- top3_preds %>% ungroup()

#########################################################

# raw top3 and iBAQ wo calibration:
save( ab_proxies , file = "data/proteomics/ALE_Proteomics_1_26_2019_ab_proxies.rdata" , compress = TRUE)
# calibrated abudances based on UPS2:
save( top3_preds , file = "data/proteomics/ALE_Proteomics_1_26_2019_top3_preds_wo_ups.rdata" , compress = TRUE)

